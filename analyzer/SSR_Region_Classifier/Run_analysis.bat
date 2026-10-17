@echo off
chcp 65001 >nul
title SSR Region Classifier

REM Switch to the directory containing this script.
cd /d "%~dp0"

REM Prefer the py launcher if available, otherwise fall back to python on PATH.
where py >nul 2>nul
if %errorlevel%==0 (
    start "" py "%~dp0ssr_region_classifier.py"
) else (
    start "" python "%~dp0ssr_region_classifier.py"
)