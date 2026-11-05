#!/usr/bin/env bash
# Gate step (round 56, task #76): prove the round witness really bracketed the round.
#
# `round_witness.sh` is tested on its own behaviour, but nothing yet checks that a
# round actually used it -- which is exactly the shape of hole this project keeps
# closing: a check that no one requires can quietly stop happening. This step asserts
# three things and is written so it can be dropped into a gate script verbatim:
#
#   1. the round has a start snapshot, and it predates the round's first gated
#      artefact (a witness started after the work would list nothing interesting);
#   2. the stop report exists and names both snapshot times, so a header-only file
#      cannot pass as "the witness ran";
#   3. coverage: every workspace path this round's commands.log says was written or
#      removed also appears in the round list. The witness is supposed to be the
#      superset of the per-command windows -- if a command-named path is missing, the
#      witness is lying, and that must be red here rather than a mystery next round.
#
# There is no built-in self-check here: the positive case is exercised by running it
# against the real round at close-out, and the refusal paths by pointing it at a round
# number that was never bracketed. If it is added to a gate script, that script's own
# mute machine should be the one that proves a skipped witness can be noticed.
#
#   usage: bash tools/gate_round_witness.sh <NN>
set -uo pipefail

root="${NEOHGT_ROOT:?source tools/env.sh first}"
tmp="${NEOHGT_TMP:?}"
cache="${NEOHGT_CACHE:?}"
round="${1:-}"
[ -n "$round" ] || { echo "usage: bash tools/gate_round_witness.sh <NN> [--selfcheck]" >&2; exit 2; }

fail=0
say() { printf '%s\n' "$*"; }
bad() { say "WITNESS GATE RED: $*"; fail=1; }

start="$tmp/.rw_${round}_before"
stop="$tmp/.rw_${round}_after"
list="$cache/logs/round${round}_files.txt"
clog="$cache/logs/commands.log"
offset_file="$tmp/.rw_${round}_clogoffset"
offset=0
[ -e "$offset_file" ] && offset=$(cat "$offset_file")
# only the bytes written after start belong to this round
clog_slice="$tmp/.gate_witness_clog_slice"
tail -c +"$((offset + 1))" "$clog" > "$clog_slice"

[ -e "$start" ] || bad "no start snapshot for round $round ($start): the round was never bracketed, so its change list is memory, not measurement"
[ -e "$list" ]  || bad "$list is missing: \`stop $round\` never ran (or refused)"
if [ -e "$list" ]; then
    # a file with only the header proves nothing ran, so both timestamps must be there
    for tag in 'start ' 'stop '; do
        grep -aq "$tag" "$list" || bad "the round list does not record a $tag time -- a header without both snapshots is not evidence"
    done
    grep -aqE 'files before=[0-9]+ after=[0-9]+' "$list" \
        || bad "the round list does not state the two snapshot sizes"
fi
if [ -e "$start" ] && [ -e "$tmp/pytest_r${round:0:2}.log" ]; then
    [ "$start" -ot "$tmp/pytest_r${round:0:2}.log" ] \
        || bad "the start snapshot is newer than the closing pass's own suite log -- the witness was started too late to cover the round"
fi

# coverage: what the per-command windows claimed vs what the round list holds.
# commands.log writes '\t\tCATEGORY\tpath', so the category is not $2 -- a first cut
# keyed on $2 and therefore reported "everything is covered" for a round whose list
# was missing entirely, which is the precise false green this step exists to catch.
missing=""
if [ -e "$list" ]; then
    missing=$(awk -F'\t' -v list="$list" '
        BEGIN { while ((getline l < list) > 0) {
                  if (l ~ /^[[:space:]]*(CREATED|MODIFIED|GONE|MOVED)[[:space:]]/) {
                      sub(/^[[:space:]]*[A-Z]+[[:space:]]+/, "", l)
                      n = split(l, a, " -> "); for (k = 1; k <= n; k++) seen[a[k]] = 1 } } }
        { cat = ""; p = ""
          for (i = 1; i < NF; i++) if ($i == "ADD" || $i == "MOD" || $i == "OUT" || $i == "GONE") { cat = $i; p = $(i+1) }
          if (cat == "") next
          sub(/ .*$/, "", p)
          if (p ~ /^neoHGT2GTDB-/ && !(p in seen) && p !~ /\.rw_/) print p }' "$clog_slice" | sort -u | head -12)
    if [ -n "$missing" ]; then
        bad "commands.log says these paths were written this round but the round list does not contain them:"
        printf '%s\n' "$missing" | sed 's/^/    /'
    else
        say "coverage: every workspace path the per-command windows named is in the round list"
    fi
fi

say "witness gate: round $round $( [ $fail = 0 ] && echo OK || echo FAILED )"
[ "$fail" = 0 ]
