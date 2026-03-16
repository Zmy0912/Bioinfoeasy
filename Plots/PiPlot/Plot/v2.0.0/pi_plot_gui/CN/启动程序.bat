@echo off
chcp 65001 > nul
echo ========================================
echo   基因动态性滑动窗口绘图程序 v2.0
echo   图形化界面版本
echo ========================================
echo.

REM 检查 Python 是否安装
python --version > nul 2>&1
if errorlevel 1 (
    echo [错误] 未检测到 Python，请先安装 Python 3.6 或更高版本
    echo.
    echo 下载地址: https://www.python.org/downloads/
    echo.
    pause
    exit /b 1
)

echo [信息] 检测到 Python 已安装
echo.

REM 检查依赖是否安装
echo [信息] 检查依赖包...
python -c "import pandas, matplotlib, numpy" > nul 2>&1
if errorlevel 1 (
    echo [警告] 检测到缺少依赖包，正在自动安装...
    echo.
    pip install -r requirements.txt
    if errorlevel 1 (
        echo.
        echo [错误] 依赖安装失败，请手动运行以下命令：
        echo         pip install pandas matplotlib numpy
        echo.
        pause
        exit /b 1
    )
    echo [信息] 依赖安装完成
    echo.
) else (
    echo [信息] 依赖包已安装
    echo.
)

REM 启动程序
echo [信息] 正在启动图形化界面...
echo.
python pi_plot_gui.py

if errorlevel 1 (
    echo.
    echo [错误] 程序运行出错
    pause
)
