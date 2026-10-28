@echo off

echo ============================================================
echo AK1 / AP5
echo 10 independent AutoDock Vina seeds
echo Same protocol as ISO
echo ============================================================

for %%S in (1001 2002 3003 4004 5005 6006 7007 8008 9009 10101) do (

    echo.
    echo ============================================================
    echo Running AP5 seed %%S
    echo ============================================================

    "C:\Users\akidi\Desktop\Docking tools\vina_win.exe" ^
    --receptor "..\receptors\AK1.pdbqt" ^
    --ligand "..\ligand\AP5.pdbqt" ^
    --center_x 25.044 ^
    --center_y 44.571 ^
    --center_z -3.708 ^
    --size_x 23.687 ^
    --size_y 23.364 ^
    --size_z 30.062 ^
    --exhaustiveness 16 ^
    --num_modes 10 ^
    --seed %%S ^
    --out "..\docking\AK1_AP5_seed%%S.pdbqt"
)

echo.
echo ============================================================
echo ALL 10 AP5 RUNS FINISHED
echo ============================================================