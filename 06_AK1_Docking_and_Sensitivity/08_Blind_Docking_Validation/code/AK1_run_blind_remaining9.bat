@echo off
setlocal

REM ============================================================
REM AK1 WHOLE-PROTEIN BLIND DOCKING
REM Remaining 9 seeds
REM Seed 1001 already completed as pilot
REM ============================================================

set "VINA=C:\Users\akidi\Desktop\Docking tools\vina_win.exe"

set "RECEPTOR=..\receptors\AK1.pdbqt"
set "LIGAND=..\ligand\isorhamnetin.pdbqt"

set "CONFIG=..\results\AK1_blind_docking\AK1_blind_box.txt"

set "EXHAUST=128"
set "MODES=20"
set "ERANGE=10"


echo.
echo ============================================================
echo AK1 BLIND DOCKING - REMAINING 9 SEEDS
echo ============================================================
echo.
echo Receptor      : %RECEPTOR%
echo Ligand        : %LIGAND%
echo Config        : %CONFIG%
echo Exhaustiveness: %EXHAUST%
echo Modes         : %MODES%
echo Energy range  : %ERANGE%
echo.
echo ============================================================


for %%S in (
2002
3003
4004
5005
6006
7007
8008
9009
10101
) do (

    echo.
    echo ============================================================
    echo BLIND DOCKING SEED %%S
    echo ============================================================
    echo.

    "%VINA%" ^
    --receptor "%RECEPTOR%" ^
    --ligand "%LIGAND%" ^
    --config "%CONFIG%" ^
    --exhaustiveness %EXHAUST% ^
    --num_modes %MODES% ^
    --energy_range %ERANGE% ^
    --seed %%S ^
    --out "..\docking\AK1_ISO_BLIND_seed%%S.pdbqt"

    if errorlevel 1 (

        echo.
        echo ERROR: seed %%S failed.
        echo Stopping batch.
        pause
        exit /b 1
    )

    echo.
    echo Seed %%S completed successfully.
)


echo.
echo ============================================================
echo ALL REMAINING 9 BLIND DOCKING RUNS COMPLETE
echo ============================================================
echo.

endlocal
pause