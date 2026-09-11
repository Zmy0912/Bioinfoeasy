@echo off
chcp 65001 > nul
echo ========================================
echo   Pi Plot - Sliding Window Analysis
echo   Graphical User Interface v2.0
echo ========================================
echo.

REM Check if Python is installed
python --version > nul 2>&1
if errorlevel 1 (
    echo [Error] Python not detected, please install Python 3.6 or higher first
    echo.
    echo Download from: https://www.python.org/downloads/
    echo.
    pause
    exit /b 1
)

echo [Info] Python detected as installed
echo.

REM Check if dependencies are installed
echo [Info] Checking dependencies...
python -c "import pandas, matplotlib, numpy" > nul 2>&1
if errorlevel 1 (
    echo [Warning] Missing dependencies detected, installing automatically...
    echo.
    pip install -r requirements.txt
    if errorlevel 1 (
        echo.
        echo [Error] Dependency installation failed, please manually run:
        echo         pip install pandas matplotlib numpy
        echo.
        pause
        exit /b 1
    )
    echo [Info] Dependencies installed successfully
    echo.
) else (
    echo [Info] Dependencies already installed
    echo.
)

REM Start the program
echo [Info] Starting graphical interface...
echo.
python pi_plot_gui.py

if errorlevel 1 (
    echo.
    echo [Error] Program execution failed
    pause
)
