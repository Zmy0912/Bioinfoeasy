@echo off
chcp 65001 >nul
echo ========================================
echo    GFF3 to BED Converter
echo ========================================
echo.
echo Starting program...
echo.

python gff3_to_bed_gui.py

if errorlevel 1 (
    echo.
    echo Program encountered an error!
    pause
)
