#!/usr/bin/env bash
# Round-level change witness.
#
# Why this exists: `tools/run.sh` compares an inventory taken before a command with
# one taken after it, so a write made *between* two commands is invisible to both
# windows. Round 54 measured it -- of that round's many documentation edits, the
# per-command logs named only the two files the commands themselves rewrote, and
# the round's change record had to be written by hand afterwards. Safety rule 16
# asks for a list of changed files per round; a list that only covers writes that
# happened to be inside a command's own window is not that list.
#
# This tool brackets a whole round: `start` snapshots the workspace, `stop`
# snapshots it again and reports every path that appeared, changed, vanished or
# (by an inferred identical-signature-and-identical-basename match) moved --
# not the write went through run.sh.
#
#   usage: bash tools/round_witness.sh start 55
#          bash tools/round_witness.sh stop  55      -> 缓存测试/logs/round55_files.txt
#
# Both snapshots are kept (nothing is deleted), and `stop` refuses without one: an
# empty witness must not be able to report "no changes" as a clean result.
set -uo pipefail

root="${NEOHGT_ROOT:?run from a shell that sourced tools/env.sh}"
tmp="${NEOHGT_TMP:?}"
cache="${NEOHGT_CACHE:?}"
action="${1:-}"
round="${2:-}"

if ! printf '%s' "$round" | grep -qE '^[1-9][0-9]{0,2}$'; then
    echo "usage: bash tools/round_witness.sh {start|stop} <round number>" >&2
    exit 2
fi
case "$action" in start|stop) ;; *)
    echo "usage: bash tools/round_witness.sh {start|stop} <round number>" >&2
    exit 2 ;;
esac

before="$tmp/.rw_${round}_before"
offset="$tmp/.rw_${round}_clogoffset"
after="$tmp/.rw_${round}_after"
out="$cache/logs/round${round}_files.txt"

# the same exclusions run.sh uses, plus the witness's own snapshots: a tool that
# reports itself as a change every time it runs teaches nobody anything
scan_excludes=(! -name '.marker.*' ! -name '.r52_inv.*' ! -name ".rw_${round}_*"
               ! -path '*/neoHGT2GTDB-临时文件/mpl/*'
               ! -path '*/xdg_cache/*' ! -path '*/__pycache__/*'
               ! -path '*/.git/*')

snapshot() {  # dest
    find "$root" -type f "${scan_excludes[@]}" -printf '%p\t%T@\t%s\n' 2>/dev/null | sort
}

if [ "$action" = start ]; then
    if [ -e "$before" ]; then
        echo "ABORT: round $round was already started at $(date -r "$before" -Is);" >&2
        echo "       restarting would hide every write made before this moment --" \
             "finish the round with \`stop\` first" >&2
        exit 3
    fi
    snapshot > "$before"
    # commands.log is append-only across rounds: without this offset the gate
    # step would compare every round's writes against this round's list and be
    # permanently red, which is how checks get switched off
    wc -c < "$cache/logs/commands.log" > "$offset"
    echo "round $round witness started: $(wc -l < "$before") files, snapshot $before," \
         "commands.log offset $(cat "$offset")"
    exit 0
fi

# stop
[ -e "$before" ] || {
    echo "ABORT: no start snapshot for round $round ($before missing)." >&2
    echo "       A witness that never started cannot report 'no changes'." >&2
    exit 4
}
snapshot > "$after"
b_count=$(wc -l < "$before"); a_count=$(wc -l < "$after")

diff_rows=$(awk -F'\t' '
    BEGIN { PROCINFO["sorted_in"] = "@ind_str_asc" }
    NR==FNR { bsig[$1] = $2 "\t" $3; next }
    { asig[$1] = $2 "\t" $3
      if (!($1 in bsig)) created[$1] = $2 "\t" $3
      else if (bsig[$1] != $2 "\t" $3) modified[$1] = 1 }
    function base(p,  q) { q = p; sub(/.*\//, "", q); return q }
    END {
      for (p in bsig) if (!(p in asig)) gone[p] = bsig[p]
      for (g in gone) for (c in created) {
          if (created[c] != gone[g]) continue
          if (base(c) != base(g)) continue
          printf "MOVED\t%s -> %s\n", g, c
          delete gone[g]; delete created[c]; break
      }
      for (c in created) printf "CREATED\t%s\n", c
      for (m in modified) printf "MODIFIED\t%s\n", m
      for (g in gone) printf "GONE\t%s\n", g
    }' "$before" "$after" | sort -t$'\t' -k1,1 -k2)

strip() { printf '%s\n' "$1" | sed "s#^$root/##"; }
{
    echo "# round $round file witness (start $(date -r "$before" -Is)," \
         "stop $(date -r "$after" -Is))"
    echo "# snapshots: $(basename "$before") ($(strip "$before")), " \
         "$(basename "$after") ($(strip "$after"))"
    echo "# workspace root: $root; files before=$b_count after=$a_count"
    [ -e "$offset" ] && echo "# commands.log offset at start: "\
        "$(cat "$offset") bytes"
    echo "# MOVED means the same basename appeared elsewhere with the same size+mtime;"
    echo "# GONE means deleted or moved outside this workspace."
    if [ -n "$diff_rows" ]; then
        printf '%s\n' "$diff_rows" | awk -F'\t' '{printf "%s\t%s\n", $1, $2}'
    else
        echo "NOTHING CHANGED -- if that is not true, the witness was started late"
    fi
} > "$out"

counts=$(printf '%s\n' "$diff_rows" | awk -F'\t' '
    NF { n[$1]++ } END { printf "created=%d modified=%d moved(inferred)=%d gone=%d",
                          n["CREATED"], n["MODIFIED"], n["MOVED"], n["GONE"] }')
counts="$counts renamed=undetected-by-design (a different basename stays CREATED+GONE; this witness claims coverage, not cleverness)"
echo "round $round witness stopped: $counts"
echo "wrote $out"
printf '%s\n' "$diff_rows" | awk -F'\t' '
    $1 == "GONE" || $1 == "MOVED" { print "  " $1 " " $2 }' | head -20
