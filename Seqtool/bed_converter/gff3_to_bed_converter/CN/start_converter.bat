@echo off
chcp 65001 >nul
echo ========================================
echo    GFF3 to BED 转换器
echo ========================================
echo.
echo 正在启动程序...
echo.

python gff3_to_bed_gui.py

if errorlevel 1 (
    echo.
    echo 程序运行出错!
    pause
)
