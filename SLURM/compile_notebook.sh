#!/usr/bin/env bash
# Nightly Lab Notebook PDF compile (WSL only).
# Replaces the former Gemini data-processing job — this script ONLY compiles the PDF.
# All placeholder processing, data analysis, notebook edits, and mirror sync are now
# handled by the Cowork/Claude nightly task (see Nightly_Cowork_Task_Prompt.md).
#
# Requirements: texlive (pdflatex) installed in WSL. No conda / py313 needed.
# Install crontab line (example: run 03:30 daily, after the Cowork task has finished editing):
#   30 3 * * * /mnt/d/Ryan\ Gustafson/PhD-Research/SLURM/compile_notebook.sh
# Note: WSL cron requires the cron service running (sudo service cron start), or launch
# this via Windows Task Scheduler calling `wsl bash -lc '/mnt/d/.../compile_notebook.sh'`.

set -u
NB_DIR="/mnt/c/Users/ryangustafson/OneDrive - University of Nevada, Reno/UNR SOM/PhD/Sarmazdeh Lab/Notebook"
LOG="/tmp/labnb_compile.log"

cd "$NB_DIR" || { echo "$(date): notebook dir not found: $NB_DIR" >> "$LOG"; exit 1; }

echo "=== compile start $(date) ===" >> "$LOG"
# Double-pass so the table of contents and cross-references resolve.
pdflatex -interaction=nonstopmode "Lab Notebook.tex" >> "$LOG" 2>&1
pdflatex -interaction=nonstopmode "Lab Notebook.tex" >> "$LOG" 2>&1

if grep -q "Output written on" "$LOG"; then
    echo "=== compile OK $(date) ===" >> "$LOG"
else
    echo "=== compile may have failed — check $LOG ($(date)) ===" >> "$LOG"
fi
