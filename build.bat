@echo off
setlocal enabledelayedexpansion

set ROOT=%~dp0
set BUILD=%ROOT%build
set SRC=%ROOT%src

rem Every model in this list is generated, validated against the dense
rem reference, and benchmarked. The build fails on any disagreement.
set MODELS=threelink fourlink elbow stanford gantry pendulums bintree satellite fourbar rrrs slidercrank quartercar springpend halfcar hmmwv hmmwv_pitman

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

echo === Building code generator ^(src\codegen^) ===
pushd "%BUILD%"
cl /nologo /Zi /EHsc /std:c++20 /Fecodegen.exe /Fo"%BUILD%\\" "%SRC%\codegen\main.cpp"
if errorlevel 1 (
    echo Code generator build failed.
    popd
    exit /b 1
)
popd

for %%M in (%MODELS%) do (
    call :valbench %%M
    if errorlevel 1 exit /b 1
)

echo === Building legacy numeric implementation ^(src\numeric^) ===
pushd "%BUILD%"
cl /nologo /Zi /std:c17 /Fenumeric.exe /Fo"%BUILD%\\" "%SRC%\numeric\main.c"
if errorlevel 1 (
    echo Legacy numeric build failed - this path is known WIP, see README.md.
)
popd

echo === Done. Artifacts and generated sources are in .\build ===
goto :eof

rem ---------------------------------------------------------------------------
rem valbench <name>: generate dynamics_<name>.c from models\<name>.mech,
rem validate it against the dense reference, and benchmark it.
rem ---------------------------------------------------------------------------
:valbench
echo === %1: generate, validate, benchmark ===
pushd "%BUILD%"
"%BUILD%\codegen.exe" "%ROOT%models\%1.mech" dynamics_%1.c
if errorlevel 1 (
    echo Code generation failed for %1.
    popd
    exit /b 1
)
cl /nologo /Zi /std:c17 /Fevalidate_%1.exe /Fo"%BUILD%\\" "%SRC%\validate\main.c" "%BUILD%\dynamics_%1.c"
if errorlevel 1 (
    echo Validation build failed for %1.
    popd
    exit /b 1
)
"%BUILD%\validate_%1.exe" "%ROOT%models\%1.mech"
if errorlevel 1 (
    echo VALIDATION FAILED: dynamics_%1.c disagrees with the dense reference.
    popd
    exit /b 1
)
cl /nologo /O2 /std:c17 /Febench_%1.exe /Fo"%BUILD%\\" "%SRC%\bench\main.c" "%BUILD%\dynamics_%1.c"
if errorlevel 1 (
    echo Benchmark build failed for %1.
    popd
    exit /b 1
)
"%BUILD%\bench_%1.exe" "%ROOT%models\%1.mech"
popd
exit /b 0
