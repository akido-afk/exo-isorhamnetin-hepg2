@echo off

echo ============================================================
echo AK1 / ISORHAMNETIN PROTONATION-STATE SENSITIVITY
echo ============================================================
echo.
echo Existing neutral baseline:
echo   Original ISO, charge 0
echo   10 seeds already completed
echo.
echo New states:
echo   7-O monoanion, charge -1
echo   7-O / 4'-O dianion, charge -2
echo.
echo ============================================================


REM ============================================================
REM STATE 1 — 7-O MONOANION
REM ============================================================

echo.
echo ============================================================
echo STATE 1: 7-O MONOANION
echo ============================================================

for %%S in (1001 2002 3003 4004 5005 6006 7007 8008 9009 10101) do (

    echo.
    echo ------------------------------------------------------------
    echo Running 7-O monoanion seed %%S
    echo ------------------------------------------------------------

    "C:\Users\akidi\Desktop\Docking tools\vina_win.exe" ^
    --receptor "..\receptors\AK1.pdbqt" ^
    --ligand "..\ligand\ISO_protonation_states\isorhamnetin_state1_7O_minus.pdbqt" ^
    --center_x 25.044 ^
    --center_y 44.571 ^
    --center_z -3.708 ^
    --size_x 23.687 ^
    --size_y 23.364 ^
    --size_z 30.062 ^
    --exhaustiveness 16 ^
    --num_modes 10 ^
    --seed %%S ^
    --out "..\docking\AK1_ISO_7Ominus_seed%%S.pdbqt"
)


REM ============================================================
REM STATE 2 — 7-O / 4'-O DIANION
REM ============================================================

echo.
echo ============================================================
echo STATE 2: 7-O / 4'-O DIANION
echo ============================================================

for %%S in (1001 2002 3003 4004 5005 6006 7007 8008 9009 10101) do (

    echo.
    echo ------------------------------------------------------------
    echo Running dianion seed %%S
    echo ------------------------------------------------------------

    "C:\Users\akidi\Desktop\Docking tools\vina_win.exe" ^
    --receptor "..\receptors\AK1.pdbqt" ^
    --ligand "..\ligand\ISO_protonation_states\isorhamnetin_state2_7O_4primeO_dianion.pdbqt" ^
    --center_x 25.044 ^
    --center_y 44.571 ^
    --center_z -3.708 ^
    --size_x 23.687 ^
    --size_y 23.364 ^
    --size_z 30.062 ^
    --exhaustiveness 16 ^
    --num_modes 10 ^
    --seed %%S ^
    --out "..\docking\AK1_ISO_dianion_seed%%S.pdbqt"
)


echo.
echo ============================================================
echo ALL PROTONATION-STATE DOCKINGS FINISHED
echo ============================================================

pause