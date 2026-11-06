# Workspace environment for neoHGT2GTDB.
#
# Source this file from any directory:
#   . <path-to-repo>/tools/env.sh
#
# It locates the code directory from this script's own location, so the
# repository may be cloned under any directory name. Every value can be
# overridden by exporting it first or by putting the assignment in
# tools/env.local (optional, git-ignored).

_neoHGT_tools_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NEOHGT_CODE="${NEOHGT_CODE:-$(cd "$_neoHGT_tools_dir/.." && pwd)}"

# Default workspace root: the directory that contains the code directory.
# Set NEOHGT_ROOT explicitly to keep outputs somewhere else.
if [ -z "${NEOHGT_ROOT:-}" ]; then
    NEOHGT_ROOT="$(cd "$NEOHGT_CODE/.." && pwd)"
fi
export NEOHGT_ROOT

export NEOHGT_TMP="${NEOHGT_TMP:-$NEOHGT_ROOT/neoHGT2GTDB-临时文件}"
export NEOHGT_CACHE="${NEOHGT_CACHE:-$NEOHGT_ROOT/neoHGT2GTDB-缓存测试}"

# Read-only GTDB release. Put the assignment in tools/env.local as well if you
# do not want to export it in every shell.
export GTDB_DIR="${GTDB_DIR:-}"
if [ -f "$NEOHGT_CODE/tools/env.local" ]; then
    . "$NEOHGT_CODE/tools/env.local"
fi

# Interpreter and DIAMOND binary. Both default to whatever is on PATH; set
# NEOHGT_PYTHON / NEOHGT_DIAMOND (or put them in tools/env.local) when the
# scientific environment lives in a conda prefix.
if [ -z "${NEOHGT_PYTHON:-}" ]; then
    NEOHGT_PYTHON="$(command -v python3 2>/dev/null || true)"
fi
export NEOHGT_PYTHON

if [ -z "${NEOHGT_DIAMOND:-}" ]; then
    NEOHGT_DIAMOND="$(command -v diamond 2>/dev/null || true)"
fi
export NEOHGT_DIAMOND
if [ -n "$NEOHGT_DIAMOND" ]; then
    PATH="$(dirname "$NEOHGT_DIAMOND"):$PATH"
    export PATH
fi

# Optional pristine upstream tree used by the zero-regression audit.
export NEOHGT_MAIN="${NEOHGT_MAIN:-}"

export PYTHONPATH="$NEOHGT_CODE${PYTHONPATH:+:$PYTHONPATH}"

# Keep third-party caches and scratch files inside the workspace.
export MPLCONFIGDIR="${MPLCONFIGDIR:-$NEOHGT_TMP/mpl}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$NEOHGT_TMP/xdg_cache}"
export TMPDIR="${NEOHGT_TMPDIR:-$NEOHGT_TMP/build}"
export NUMEXPR_CACHE_DIR="${NUMEXPR_CACHE_DIR:-$NEOHGT_TMP/numexpr}"
export PYTHONDONTWRITEBYTECODE=1
export PIP_CACHE_DIR="${PIP_CACHE_DIR:-$NEOHGT_TMP/pip_cache}"
export PIP_DISABLE_PIP_VERSION_CHECK=1

mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME" "$TMPDIR" "$PIP_CACHE_DIR" \
         "$NEOHGT_CACHE/logs" "$NEOHGT_TMP/build" "$NEOHGT_TMP/search"

if [ -z "${NEOHGT_PYTHON:-}" ]; then
    echo "WARNING: python3 was not found on PATH; set NEOHGT_PYTHON." >&2
fi
if [ -z "${GTDB_DIR:-}" ]; then
    echo "NOTE: GTDB_DIR is empty; pass --gtdb-dir or set it in tools/env.local." >&2
fi
