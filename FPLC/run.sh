#!/usr/bin/env bash
# Launch the FPLC analyzer app using the Windows conda streamlit environment
STREAMLIT="/c/Users/ryangustafson/AppData/Local/miniconda3/envs/streamlit/Scripts/streamlit"
APP_DIR="$(cd "$(dirname "$0")" && pwd)"
"$STREAMLIT" run "$APP_DIR/app.py"
