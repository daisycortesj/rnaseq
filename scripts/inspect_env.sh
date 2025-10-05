
#!/usr/bin/env bash
set -euo pipefail

echo "[info] PATH: $PATH"
echo "[info] Conda env: ${CONDA_DEFAULT_ENV:-none}"

if command -v STAR >/dev/null 2>&1; then
  echo -n "[ok] STAR version: "
  STAR --version || true
else
  echo "[warn] STAR not on PATH"
fi

if command -v Trinity >/dev/null 2>&1; then
  echo -n "[ok] Trinity version: "
  Trinity --version || true
else
  echo "[warn] Trinity not on PATH"
fi
