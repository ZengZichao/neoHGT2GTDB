#!/usr/bin/env bash
# Audited command runner.
#
#   usage: run.sh "<intent>" "<impact on files>" -- <command> [args...]
#
# Prints cwd, the full command, the declared impact, runs the command, then
# appends an entry to $NEOHGT_CACHE/logs/commands.log and lists every file
# created, modified, renamed or removed inside $NEOHGT_ROOT as a result of that
# command -- by comparing two inventories of the tree, not by asking find which
# files are newer (a removed file cannot be newer than anything).
#
# Two pre-flight refusals: the destructive patterns banned by the safety rules
# (rule 14), and redirections whose target resolves outside $NEOHGT_ROOT.
#
# A redirection check alone does NOT cover a tool that writes somewhere else by
# its own default -- matplotlib's cache, a tempdir, a ~/.config file never
# appear in the command line. So the runner also scans $NEOHGT_WATCH (external
# locations, colon separated) after the command and reports anything that
# changed there. With NEOHGT_WATCH_STRICT=1 such a write fails the run.

set -uo pipefail

intent="$1"; shift
impact="$1"; shift
[ "${1:-}" = "--" ] && shift
if [ "$#" -eq 0 ]; then
    echo "run.sh: no command given" >&2
    exit 2
fi
cmd="$*"
# The command runs from the argument vector. It used to run through
# `eval "$*"`, which drops every quote: `-- python -c "print('x')"` reached the
# shell as `python -c print('x')` and died on a syntax error (round 52,
# measured). The joined text is still what the pre-flight checks below read --
# a ban has to look at what a human typed -- but the executed and logged form
# is the vector.
argv=("$@")
cmd_shown=$(printf '%q ' "${argv[@]}")

root=${NEOHGT_ROOT:?source tools/env.sh first}
log="$root/neoHGT2GTDB-缓存测试/logs/commands.log"

# Banned destructive patterns (safety rule 14).
if printf '%s' "$cmd" | grep -qE 'rm[[:space:]]+-[a-zA-Z]*r[a-zA-Z]*f|rm[[:space:]]+-[a-zA-Z]*f[a-zA-Z]*r|find[[:space:]].*-delete|git[[:space:]]+(clean|reset[[:space:]]+--hard|checkout[[:space:]]+--)|:\(\)[[:space:]]*\{[[:space:]]*:[[:space:]]*\|[[:space:]]*:\|[[:space:]]*\}'; then
    echo "run.sh: REFUSED — command matches a banned destructive pattern" >&2
    exit 3
fi
# Only absolute-path redirections can escape the workspace; relative ones land
# in the (contained) cwd, and arrows inside strings ("->") must not trip this.
esc=$(printf '%s' "$cmd" | grep -oE '[0-9]*>{1,2}[[:space:]]*(/[[:graph:]]+|\$HOME[[:graph:]]*|~[[:graph:]]*)' \
        | sed -E 's/^[0-9]*>{1,2}[[:space:]]*//')
for tgt in $esc; do
    case "$tgt" in
        \$HOME*|~*) echo "run.sh: REFUSED — redirection targets \$HOME ($tgt)" >&2; exit 3 ;;
    esac
    norm=$(realpath -m "$tgt" 2>/dev/null || printf '%s' "$tgt")
    case "$norm" in
        "$root"/*) ;;
        *) echo "run.sh: REFUSED — redirection target outside $root: $tgt -> $norm" >&2; exit 3 ;;
    esac
done

watch=${NEOHGT_WATCH:-"$HOME/.cache:$HOME/.config:$HOME/.matplotlib:/tmp:/var/tmp:/dev/shm"}

# The workspace inventory is compared before/after, not filtered by mtime.
# `find -newer` can only answer "which files are newer?", and a file that was
# deleted is not there to answer it; a rename keeps its mtime, so the new path
# is not newer either. Round 52 measured it (r99_runner_change_list_blindness.py,
# with a write as positive control): the change line below printed "(none)" for
# both a `mv` and a removal -- the two operations safety rule 14 needs listed.
scan_excludes=(! -name '.marker.*' ! -name '.r52_inv.*'
               ! -path '*/neoHGT2GTDB-临时文件/mpl/*'
               ! -path '*/xdg_cache/*' ! -path '*/__pycache__/*')
inventory() {
    find "$root" -type f "${scan_excludes[@]}" -printf '%p\t%T@\t%s\n' 2>/dev/null | sort
}

marker=$(mktemp "$root/neoHGT2GTDB-临时文件/.marker.XXXXXX")
inv_before=$(mktemp "$root/neoHGT2GTDB-临时文件/.r52_inv.XXXXXX")
inv_after=$(mktemp "$root/neoHGT2GTDB-临时文件/.r52_inv.XXXXXX")
echo "----------------------------------------------------------------"
echo "intent : $intent"
echo "impact : $impact"
echo "cwd    : $(pwd)"
echo "command: $cmd_shown"
echo "----------------------------------------------------------------"

inventory > "$inv_before"
start=$(date -u +%Y-%m-%dT%H:%M:%SZ)
"${argv[@]}"
rc=$?
end=$(date -u +%Y-%m-%dT%H:%M:%SZ)
inventory > "$inv_after"

external=""
for dir in $watch; do
    [ -d "$dir" ] || continue
    hit=$(find "$dir" -maxdepth 3 -type f -newer "$marker" 2>/dev/null | head -20)
    [ -n "$hit" ] && external="$external\n$hit"
done
if [ -n "$external" ]; then
    echo "== EXTERNAL writes outside $root (declared in the impact line?)"
    printf "%b\n" "$external" | sed '/^$/d'
fi

# CREATED / MODIFIED / MOVED / AMBIG / GONE, one category per line, path stripped
# of the $root prefix. MOVED is an inference: a path that vanished and one that
# appeared with the same mtime and size -- labelled as inferred because a
# signature cannot prove identity (two empty files look alike). Since round 54 the
# inference also requires the file *name* to match, and says AMBIG instead of
# guessing when several created paths share the signature: moving a whole tree
# makes every gone path signature-identical to every created one of the same size,
# and the first version printed an arbitrary permutation of that as fact (192 rows
# in one command, pairing docs/A with a different tree's docs/B).
diff_rows=$(awk -F'\t' '
    NR==FNR { bsig[$1] = $2 "\t" $3; next }
    { asig[$1] = $2 "\t" $3
      if (!($1 in bsig)) { created[$1] = $2 "\t" $3 }
      else if (bsig[$1] != $2 "\t" $3) { modified[$1] = 1 } }
    function base(p,  q) { q = p; sub(/.*\//, "", q); return q }
    END {
      for (p in bsig) if (!(p in asig)) { gone[p] = bsig[p]; n++ ; gk[n] = p }
      for (i = 1; i <= n; i++) {
          g = gk[i]
          if (!(g in gone)) continue
          cand = 0; same = 0; pick = ""
          for (c in created) {
              if (created[c] != gone[g]) continue
              cand++
              if (base(c) == base(g)) { same++; pick = c }
          }
          if (same == 1) {
              printf "MOVED\t%s -> %s\n", g, pick
              delete gone[g]; delete created[pick]; continue
          }
          if (cand == 1 && same == 0) {
              for (c in created) if (created[c] == gone[g]) {
                  printf "MOVED\t%s -> %s (renamed)\n", g, c
                  delete gone[g]; delete created[c] }
              continue
          }
          if (cand > 1) {
              printf "AMBIG\t%s (%s created paths share its size+mtime; " \
                     "not claiming which)\n", g, cand
              delete gone[g]
          }
      }
      for (c in created) printf "CREATED\t%s\n", c
      for (m in modified) printf "MODIFIED\t%s\n", m
      for (g in gone) printf "GONE\t%s (deleted, or moved outside the workspace)\n", g
    }' "$inv_before" "$inv_after" | sort -t$'\t' -k1,1 -k2)

n_created=$(printf '%s\n' "$diff_rows" | grep -c $'^CREATED\t')
n_modified=$(printf '%s\n' "$diff_rows" | grep -c $'^MODIFIED\t')
n_moved=$(printf '%s\n' "$diff_rows" | grep -c $'^MOVED\t')
n_ambig=$(printf '%s\n' "$diff_rows" | grep -c $'^AMBIG\t')
n_gone=$(printf '%s\n' "$diff_rows" | grep -c $'^GONE\t')

echo "== changed files under $root (before/after inventory, so a delete or a rename is visible too) =="
if [ -z "$diff_rows" ]; then
    echo "(none)"
else
    printf 'created=%s modified=%s moved(inferred)=%s moved(ambiguous)=%s gone=%s\n' \
        "$n_created" "$n_modified" "$n_moved" "$n_ambig" "$n_gone"
    for cat in GONE AMBIG MOVED CREATED MODIFIED; do
        rows=$(printf '%s\n' "$diff_rows" | grep "^${cat}	" | sed "s|^${cat}	||" | sed "s|$root/||g")
        [ -n "$rows" ] || continue
        printf '%s\n' "$rows" | head -20 | sed "s/^/  $cat  /"
        total=$(printf '%s\n' "$rows" | wc -l)
        [ "$total" -gt 20 ] && echo "  $cat ... ($total in total, first 20 shown)"
    done
fi

{
    printf '%s\t%s\t%s\texit=%s\tcwd=%s\t%s\n' "$start" "$end" "$intent" "$rc" "$(pwd)" "$cmd_shown"
    if [ -n "$external" ]; then
        printf '\t\tEXTERNAL\t'
        printf "%b\n" "$external" | sed '/^$/d' | paste -sd, -
    fi
    printf '%s\n' "$diff_rows" | sed "s|$root/||" \
        | sed "s|^CREATED\\t|\\t\\tOUT\\t|;s|^MODIFIED\\t|\\t\\tMOD\\t|;s|^GONE\\t|\\t\\tGONE\\t|;s|^MOVED\\t|\\t\\tMOVED\\t|;s|^AMBIG\\t|\\t\\tAMBIG\\t|"
    echo ""
} >> "$log"

rm -f "$marker" "$inv_before" "$inv_after"
if [ -n "$external" ] && [ "${NEOHGT_WATCH_STRICT:-0}" = "1" ]; then
    echo "run.sh: FAILED — the command wrote outside $root while " >&2
    echo "        NEOHGT_WATCH_STRICT=1 (see the EXTERNAL list above)" >&2
    exit 4
fi
exit $rc
