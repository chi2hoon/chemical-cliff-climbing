#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_DIR"

PYTHON_BIN="${PYTHON_BIN:-python3}"

if [ ! -d ".venv" ]; then
  echo "[setup] Creating virtual environment (.venv)"
  "$PYTHON_BIN" -m venv .venv
fi

echo "[setup] Activating virtual environment"
# shellcheck disable=SC1091
source .venv/bin/activate

echo "[setup] Upgrading pip/setuptools/wheel"
python -m pip install --upgrade pip setuptools wheel

echo "[setup] Installing requirements"
pip install -r requirements.txt

echo "[run] Launching Streamlit app"
exec streamlit run app.py


