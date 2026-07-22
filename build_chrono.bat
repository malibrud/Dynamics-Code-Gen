@echo off
setlocal enabledelayedexpansion

rem Builds and runs the Chrono cross-check (src\chrono): every .mech model is
rem rebuilt as a Project Chrono multibody system and its accelerations are
rem compared against the dense RBDA reference. Separate from build.bat so the
rem main build has no Chrono dependency.
rem
rem Requires the Chrono SDK and Eigen; override the default locations with
rem the CHRONO_DIR and EIGEN_DIR environment variables.

set ROOT=%~dp0
set BUILD=%ROOT%build
set SRC=%ROOT%src

if not defined CHRONO_DIR set CHRONO_DIR=C:\dev\chrono
if not defined EIGEN_DIR set EIGEN_DIR=C:\dev\eigen-5.0.0

if not exist "%CHRONO_DIR%\include\chrono\physics\ChSystemNSC.h" (
    echo Chrono SDK not found at %CHRONO_DIR% - set CHRONO_DIR.
    exit /b 1
)
if not exist "%EIGEN_DIR%\Eigen" (
    echo Eigen not found at %EIGEN_DIR% - set EIGEN_DIR.
    exit /b 1
)

set MODELS=threelink fourlink elbow stanford gantry pendulums bintree satellite fourbar rrrs slidercrank quartercar springpend halfcar hmmwv hmmwv_pitman
rem Optional: pass model names as arguments to check just those.
if not "%~1"=="" set MODELS=%*

if not exist "%BUILD%" mkdir "%BUILD%"

rem If cl.exe isn't already on PATH (e.g. not running from a Developer
rem Command Prompt), locate Visual Studio via vswhere and init the x64 tools.
where cl.exe >nul 2>nul
if errorlevel 1 (
    set "PF86=%ProgramFiles(x86)%"
    set "VSWHERE=!PF86!\Microsoft Visual Studio\Installer\vswhere.exe"
    if not exist "!VSWHERE!" (
        echo Could not find vswhere.exe. Install Visual Studio with the C++
        echo workload, or run this from a Developer Command Prompt.
        exit /b 1
    )
    for /f "usebackq tokens=*" %%i in (`"!VSWHERE!" -latest -property installationPath`) do set "VSPATH=%%i"
    if not defined VSPATH (
        echo Could not locate a Visual Studio installation via vswhere.
        exit /b 1
    )
    call "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat" >nul
)

rem The checker links each model's generated dynamics, so it needs the code
rem generator. Built here rather than assuming build.bat has run, to keep the
rem two scripts independent.
echo === Building code generator ^(src\codegen^) ===
pushd "%BUILD%"
cl /nologo /O2 /EHsc /std:c++20 /Fecodegen.exe /Fo"%BUILD%\\" "%SRC%\codegen\main.cpp"
if errorlevel 1 (
    echo Code generator build failed.
    popd
    exit /b 1
)
popd

rem Chrono_core.dll etc. live in the SDK's bin directory.
set PATH=%CHRONO_DIR%\bin;%PATH%

for %%M in (%MODELS%) do (
    call :check %%M
    if errorlevel 1 exit /b 1
)

echo === Done. All models agree with the Chrono engine. ===
goto :eof

rem ---------------------------------------------------------------------------
rem check <name>: generate dynamics_<name>.c, build a chrono checker linking
rem it, and run the cross-check + simulation comparison against Chrono.
rem ---------------------------------------------------------------------------
:check
echo === %1: Chrono cross-check ===
pushd "%BUILD%"
"%BUILD%\codegen.exe" "%ROOT%models\%1.mech" dynamics_%1.c
if errorlevel 1 (
    echo Code generation failed for %1.
    popd
    exit /b 1
)
rem Per-model object directory: main.obj is compiled once per model here, and
rem would otherwise collide with build.bat's own main.obj in %BUILD%.
if not exist "%BUILD%\obj_chrono_%1" mkdir "%BUILD%\obj_chrono_%1"
rem /arch:AVX2 and /Zc:__cplusplus are required by the Chrono SDK build
rem (see INTERFACE_COMPILE_OPTIONS in ChronoTargets.cmake): Chrono was
rem compiled with AVX2, so Eigen objects crossing the DLL boundary use
rem 32-byte alignment -- mixing architectures corrupts the heap.
cl /nologo /O2 /MD /EHsc /std:c++17 /bigobj /arch:AVX2 /Zc:__cplusplus ^
   /DNOMINMAX /DWIN32_LEAN_AND_MEAN /D_CRT_SECURE_NO_WARNINGS ^
   /I"%CHRONO_DIR%\include" /I"%CHRONO_DIR%\include\chrono_thirdparty\yaml-cpp\include" /I"%EIGEN_DIR%" ^
   /Fechrono_check_%1.exe /Fo"%BUILD%\obj_chrono_%1\\" "%SRC%\chrono\main.cpp" "%BUILD%\dynamics_%1.c" ^
   /link /LIBPATH:"%CHRONO_DIR%\lib" Chrono_core.lib yaml-cpp.lib
if errorlevel 1 (
    echo Chrono cross-check build failed for %1.
    popd
    exit /b 1
)
popd
"%BUILD%\chrono_check_%1.exe" "%ROOT%models\%1.mech"
rem 'neq 0' rather than 'errorlevel 1': a crashed process returns a negative
rem code, which 'if errorlevel 1' would treat as success.
if !errorlevel! neq 0 (
    echo CHRONO CROSS-CHECK FAILED for %1.
    exit /b 1
)
exit /b 0
