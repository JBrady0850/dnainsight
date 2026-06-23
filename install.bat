@echo off
setlocal EnableDelayedExpansion

echo.
echo ============================================================
echo   DNAInsight Installer - Windows
echo ============================================================
echo.

:: ============================================================
:: STEP 0: Check / Auto-install Python
:: ============================================================

python --version >nul 2>&1
if not errorlevel 1 goto :python_ok

echo [WARN] Python not found in PATH.
echo [AUTO] Attempting to install Python 3.12 via Windows Package Manager (winget)...
echo.
winget install --id Python.Python.3.12 --silent --accept-package-agreements --accept-source-agreements >nul 2>&1
if not errorlevel 1 (
    echo [OK] Python installed via winget.
    goto :check_path_reload
)

echo [AUTO] winget install failed or not available. Downloading Python installer...
powershell -NoProfile -Command ^
    "[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; " ^
    "try { " ^
    "  (New-Object Net.WebClient).DownloadFile(" ^
    "    'https://www.python.org/ftp/python/3.12.7/python-3.12.7-amd64.exe'," ^
    "    '%TEMP%\py_setup.exe'" ^
    "  ); exit 0 " ^
    "} catch { exit 1 }" >nul 2>&1
if errorlevel 1 (
    echo.
    echo [ERROR] Could not auto-download Python.
    echo.
    echo Please install Python 3.10+ manually:
    echo   https://www.python.org/downloads/
    echo   During install: check "Add Python to PATH"
    echo.
    echo Then re-run this installer.
    pause
    exit /b 1
)
echo [AUTO] Running Python installer silently...
%TEMP%\py_setup.exe /quiet InstallAllUsers=0 PrependPath=1 Include_test=0 Include_pip=1
del /f /q %TEMP%\py_setup.exe >nul 2>&1

:check_path_reload
:: Try to pick up Python from known install paths without requiring a new shell
set PYTHON_FOUND=0
for %%V in (313 312 311 310) do (
    if exist "%LOCALAPPDATA%\Programs\Python\Python%%V\python.exe" (
        set "PATH=%LOCALAPPDATA%\Programs\Python\Python%%V;%LOCALAPPDATA%\Programs\Python\Python%%V\Scripts;%PATH%"
        set PYTHON_FOUND=1
        goto :python_ok
    )
    if exist "C:\Program Files\Python%%V\python.exe" (
        set "PATH=C:\Program Files\Python%%V;C:\Program Files\Python%%V\Scripts;%PATH%"
        set PYTHON_FOUND=1
        goto :python_ok
    )
)

:: If still not found, instruct user to reopen
echo.
echo [INFO] Python has been installed but PATH has not refreshed in this window.
echo.
echo Please CLOSE this window and run install.bat again to complete setup.
echo.
pause
exit /b 0

:python_ok
for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYVER=%%i
echo [OK] Python %PYVER% detected.

:: ============================================================
:: STEP 1: Install / upgrade pip
:: ============================================================
echo.
echo [1/4] Checking pip...
python -m pip --version >nul 2>&1
if errorlevel 1 (
    echo [AUTO] pip not found. Installing via ensurepip...
    python -m ensurepip --upgrade >nul 2>&1
)
python -m pip install --upgrade pip --quiet >nul 2>&1
echo [OK] pip ready.

:: ============================================================
:: STEP 2: Install Python dependencies
:: ============================================================
echo.
echo [2/4] Installing Python dependencies...
pip install -r requirements.txt --quiet
if errorlevel 1 (
    echo [ERROR] Dependency installation failed. Check your internet connection.
    pause
    exit /b 1
)
echo [OK] Dependencies installed.

:: ============================================================
:: STEP 3: Build bundled SNP reference
:: ============================================================
echo.
echo [3/4] Building bundled SNP reference database...
python data\build_reference.py
if errorlevel 1 (
    echo [ERROR] Failed to build SNP reference.
    pause
    exit /b 1
)
echo [OK] SNP reference ready.

:: ============================================================
:: STEP 4: Create launcher and directories
:: ============================================================
echo.
echo [4/4] Creating launch script and directories...

echo @echo off > launch.bat
echo cd /d "%%~dp0" >> launch.bat
echo python app.py %%* >> launch.bat

if not exist uploads mkdir uploads
if not exist reports_output mkdir reports_output
if not exist db mkdir db
if not exist data mkdir data

echo [OK] Setup complete.

echo.
echo ============================================================
echo   Installation Complete!
echo ============================================================
echo.
echo To launch DNAInsight:
echo   Double-click launch.bat
echo   -- or -- python app.py
echo.
echo DNAInsight opens in your browser at http://127.0.0.1:5050
echo.
echo TIP: Open DNAInsight and use Settings > Database to update
echo      ClinVar annotations monthly for the best results.
echo.

set /p LAUNCH="Launch DNAInsight now? (Y/N): "
if /i "!LAUNCH!"=="Y" (
    start "" python app.py
)

pause
