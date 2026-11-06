"""Test-session environment bootstrap.

CI runners (and any checkout without tools/env.sh) used to fail before the
first assertion: tests index os.environ['NEOHGT_CACHE'] directly, and the
workspace guard compared pytest's /tmp output against NEOHGT_ROOT. Defaulting
the four variables into one scratch workspace keeps a bare checkout green or
honestly skipped, without masking anything env.sh sets.
"""

import os
import tempfile


def _bootstrap():
    root = os.environ.get('NEOHGT_ROOT')
    if not root:
        # a disposable workspace: the guard can then contain writes, and the
        # read-only-input rule still applies on top of it
        root = tempfile.mkdtemp(prefix='neohgt_test_ws_')
        os.environ['NEOHGT_ROOT'] = root
    code = os.environ.get('NEOHGT_CODE')
    if not code:
        # conftest.py lives at <root>/neoHGT/tests/: the code root is three
        # levels up (the previous two left NEOHGT_CODE at the package dir,
        # which pointed every ledger test at neoHGT/docs/...)
        os.environ['NEOHGT_CODE'] = os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__))))
    for var, sub in (('NEOHGT_CACHE', 'cache'), ('NEOHGT_TMP', 'tmp')):
        if not os.environ.get(var):
            os.environ[var] = os.path.join(root, sub)
        # the helper and the guard both expect these to exist; creating them
        # here keeps the scratch tree inside the workspace even when the host
        # already exports TMPDIR
        os.makedirs(os.environ[var], exist_ok=True)
    if not os.environ.get('TMPDIR'):
        # keep pytest temp dirs inside the workspace so guard_write sees them
        tmp = os.path.join(os.environ['NEOHGT_TMP'], 'build')
        os.makedirs(tmp, exist_ok=True)
        os.environ['TMPDIR'] = tmp
    # byte-code caches inside a source tree made one gate mistake an
    # environment artifact for a product regression
    os.environ.setdefault('PYTHONDONTWRITEBYTECODE', '1')


_bootstrap()
