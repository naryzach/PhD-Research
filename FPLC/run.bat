@echo off
REM Launch the FPLC Analyzer using the streamlit conda environment
set CONDA_ROOT=C:\Users\ryangustafson\AppData\Local\miniconda3
set STREAMLIT=%CONDA_ROOT%\envs\streamlit\Scripts\streamlit.exe
set APP=\\wsl.localhost\ubuntu\home\ryangustafson\Documents\GitHubProj\PhD-Research\FPLC\app.py

echo Starting FPLC Analyzer...
"%STREAMLIT%" run "%APP%"
pause
