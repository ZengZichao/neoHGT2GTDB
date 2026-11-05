# Workspace environment for the neoHGT2GTDB adaptation.
# Source it from anywhere -- it locates the workspace relative to itself:
#   . <workspace>/neoHGT2GTDB-项目代码/tools/env.sh
#
# Every writable path below lives inside NEOHGT_ROOT; third-party caches are
# redirected into NEOHGT_TMP so that nothing is written to /tmp or $HOME.

# workspace root = the directory that contains this project folder
export NEOHGT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export NEOHGT_CODE=$NEOHGT_ROOT/neoHGT2GTDB-项目代码
export NEOHGT_TMP=$NEOHGT_ROOT/neoHGT2GTDB-临时文件
export NEOHGT_CACHE=$NEOHGT_ROOT/neoHGT2GTDB-缓存测试

# Read-only input dataset (never written to; enforced by neoHGT.gtdb guards).
# Point GTDB_DIR at your local GTDB release copy before sourcing, or put the
# assignment (a path, never a secret) in tools/env.local -- it is sourced
# below and git-ignored.
export GTDB_DIR="${GTDB_DIR:-}"
if [ -f "$NEOHGT_CODE/tools/env.local" ]; then
    . "$NEOHGT_CODE/tools/env.local"
fi

# Interpreter and DIAMOND binary come from pre-existing conda installations;
# nothing is downloaded or installed by this project.
export NEOHGT_PYTHON=${NEOHGT_PYTHON:-/opt/anaconda3/bin/python}
# the pristine upstream tree the reachability audit diffs against
export NEOHGT_MAIN=${NEOHGT_MAIN:-$NEOHGT_ROOT/neoHGT-main}

export NEOHGT_DIAMOND=${NEOHGT_DIAMOND:-/opt/anaconda3/envs/diamond-2.1.10/bin/diamond}
export PATH="$(dirname "$NEOHGT_DIAMOND"):$PATH"

export PYTHONPATH=$NEOHGT_CODE${PYTHONPATH:+:$PYTHONPATH}

# Cache/temp redirection (constraint: no writes outside NEOHGT_ROOT).
export MPLCONFIGDIR=$NEOHGT_TMP/mpl
export XDG_CACHE_HOME=$NEOHGT_TMP/xdg_cache
export TMPDIR=$NEOHGT_TMP/build
export NUMEXPR_CACHE_DIR=$NEOHGT_TMP/numexpr
export PYTHONDONTWRITEBYTECODE=1
# Round 52: IT-7 attributes a fresh file under a cache namespace this stack uses
# to this stack. pip's "is there a newer pip" self-check lives in ~/.cache/pip, so
# redirect pip's whole cache in here and switch the version check off: nothing we
# launch can then write there, which is what makes the declared-observer exemption
# in verify_workspace.py a measurement rather than an excuse.
export PIP_CACHE_DIR=$NEOHGT_TMP/pip_cache
export PIP_DISABLE_PIP_VERSION_CHECK=1

mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$TMPDIR" "$PIP_CACHE_DIR" \
         "$NEOHGT_CACHE/logs" "$NEOHGT_TMP/build" "$NEOHGT_TMP/search"
