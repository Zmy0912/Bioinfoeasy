@echo off
chcp 65001 >nul
title SSR 区域分类分析工具

REM 切换到脚本所在目录
cd /d "%~dp0"

REM 尝试用 python（py launcher 优先，否则回退到 PATH 中的 python）
where py >nul 2>nul
if %errorlevel%==0 (
    start "" py "%~dp0ssr_region_classifier.py"
) else (
    start "" python "%~dp0ssr_region_classifier.py"
)