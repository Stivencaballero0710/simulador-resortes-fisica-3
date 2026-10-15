@echo off
cd /d "%~dp0"
where py >nul 2>&1
if not errorlevel 1 (
    py -3 -m venv .venv
) else (
    python -m venv .venv
)
if errorlevel 1 goto fallo
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto fallo
echo.
echo Instalacion terminada. Ejecute verificar.bat y luego abrir_simulador.bat.
pause
exit /b 0
:fallo
echo.
echo La instalacion fallo. Revise Python 3.11 o posterior y la conexion a internet.
pause
exit /b 1

