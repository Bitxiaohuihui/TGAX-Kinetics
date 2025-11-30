@echo off
REM =================================================================
REM  TGAX Kinetics PyInstaller Build Script
REM =================================================================
REM  This script automates the process of building the TGAX Kinetics.exe
REM  It must be placed in the same folder as:
REM  - TGAX_Kinetics.py
REM  - splash.png
REM  - BIT_Kinetics_Icon_Tight.ico
REM =================================================================

echo.
echo [INFO] Starting the TGAX Kinetics application build process...
echo [INFO] This may take a few minutes. Please wait.
echo.

REM --- Check for Python ---
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python is not found in your system's PATH.
    echo Please install Python and make sure it's added to the PATH.
    goto :error
)

REM --- Check for PyInstaller ---
pip show pyinstaller >nul 2>&1
if %errorlevel% neq 0 (
    echo [INFO] PyInstaller not found. Attempting to install it now...
    pip install pyinstaller
    if %errorlevel% neq 0 (
        echo [ERROR] Failed to install PyInstaller. Please install it manually using 'pip install pyinstaller'.
        goto :error
    )
    echo [SUCCESS] PyInstaller has been installed.
)

REM --- Run the PyInstaller command ---
echo.
echo [INFO] Running PyInstaller to build the .exe file...
pyinstaller --name "TGAX Kinetics" --onefile --windowed --icon="BIT_Kinetics_Icon_Tight.ico" --add-data "splash.png;." --add-data "BIT_Kinetics_Icon_Tight.ico;." TGAX_Kinetics.py

REM --- Check build result ---
if %errorlevel% neq 0 (
    echo [ERROR] PyInstaller failed to build the application. Please check the messages above for errors.
    goto :error
)

echo.
echo =================================================================
echo [SUCCESS] Build completed successfully!
echo.
echo Your application "TGAX Kinetics.exe" can be found in the 'dist' folder.
echo =================================================================
echo.
goto :end

:error
echo.
echo Build failed. Please correct the errors and try again.
pause

:end
pause