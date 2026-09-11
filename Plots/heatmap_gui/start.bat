@echo off
cd /d "%~dp0"
python heatmap_gui_v3.py
if errorlevel 1 pause
