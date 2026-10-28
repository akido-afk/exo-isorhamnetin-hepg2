@echo off
setlocal

REM ============================================================
REM AK1 RECEPTOR-CONFORMATION SENSITIVITY
REM Neutral isorhamnetin against three selected 7X7S conformers
REM ============================================================

set "VINA=C:\Users\akidi\Desktop\Docking tools\vina_win.exe"
set "LIGAND=..\ligand\isorhamnetin.pdbqt"

REM ============================================================
REM SAME 1Z83 NUCLEOTIDE-REGION SEARCH BOX
REM ============================================================

set "CX=25.044"
set "CY=44.571"
set "CZ=-3.708"

set "SX=23.687"
set "SY=23.364"
set "SZ=30.062"

set "EXHAUST=16"
set "MODES=10"

echo.
echo ============================================================
echo AK1 7X7S ENSEMBLE DOCKING
echo ============================================================
echo Ligand: neutral isorhamnetin
echo Exhaustiveness: %EXHAUST%
echo Modes: %MODES%
echo Seeds: 10 per receptor
echo Total docking runs: 30
echo ============================================================


REM ============================================================
REM MODEL 20 - LEAST OPEN
REM ============================================================

echo.
echo ============================================================
echo MODEL 20 - LEAST OPEN
echo ============================================================

for %%S in (1001 2002 3003 4004 5005 6006 7007 8008 9009 10101) do (

    echo.
    echo ------------------------------------------------------------
    echo Model 20 - seed %%S
    echo ------------------------------------------------------------

    "%VINA%" ^
    --receptor "..\receptors\AK1_7X7S_model20_least_open.pdbqt" ^
    --ligand "%LIGAND%" ^
    --center_x %CX% ^
    --center_y %CY% ^
    --center_z %CZ% ^
    --size_x %SX% ^
    --size_y %SY% ^
    --size_z %SZ% ^
    --exhaustiveness %EXHAUST% ^
    --num_modes %MODES% ^
    --seed %%S ^
    --out "..\docking\AK1_7X7S_model20_least_open_ISO_seed%%S.pdbqt"
)


REM ============================================================
REM MODEL 8 - ENSEMBLE MEDOID
REM ============================================================

echo.
echo ============================================================
echo MODEL 8 - ENSEMBLE MEDOID
echo ============================================================

for %%S in (1001 2002 3003 4004 5005 6006 7007 8008 9009 10101) do (

    echo.
    echo ------------------------------------------------------------
    echo Model 8 - seed %%S
    echo ------------------------------------------------------------

    "%VINA%" ^
    --receptor "..\receptors\AK1_7X7S_model8_ensemble_medoid.pdbqt" ^
    --ligand "%LIGAND%" ^
    --center_x %CX% ^
    --center_y %CY% ^
    --center_z %CZ% ^
    --size_x %SX% ^
    --size_y %SY% ^
    --size_z %SZ% ^
    --exhaustiveness %EXHAUST% ^
    --num_modes %MODES% ^
    --seed %%S ^
    --out "..\docking\AK1_7X7S_model8_ensemble_medoid_ISO_seed%%S.pdbqt"
)


REM ============================================================
REM MODEL 1 - MOST OPEN
REM ============================================================

echo.
echo ============================================================
echo MODEL 1 - MOST OPEN
echo ============================================================

for %%S in (1001 2002 3003 4004 5005 6006 7007 8008 9009 10101) do (

    echo.
    echo ------------------------------------------------------------
    echo Model 1 - seed %%S
    echo ------------------------------------------------------------

    "%VINA%" ^
    --receptor "..\receptors\AK1_7X7S_model1_most_open.pdbqt" ^
    --ligand "%LIGAND%" ^
    --center_x %CX% ^
    --center_y %CY% ^
    --center_z %CZ% ^
    --size_x %SX% ^
    --size_y %SY% ^
    --size_z %SZ% ^
    --exhaustiveness %EXHAUST% ^
    --num_modes %MODES% ^
    --seed %%S ^
    --out "..\docking\AK1_7X7S_model1_most_open_ISO_seed%%S.pdbqt"
)


echo.
echo ============================================================
echo ALL 30 ENSEMBLE DOCKING RUNS FINISHED
echo ============================================================
echo.

endlocal
pause