@echo off
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
    echo Ejecute primero instalar.bat.
    pause
    exit /b 1
)
".venv\Scripts\python.exe" simulador.py --interactivo
set "RESULTADO=%ERRORLEVEL%"
echo.
if not "%RESULTADO%"=="0" echo La ejecucion termino con un error. Revise el mensaje anterior.
pause
exit /b %RESULTADO%

