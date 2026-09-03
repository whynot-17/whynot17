# STAT3 / DINP MD protocol template

Selected RCSB entry: 6NJS
Requested production length: 100 ns
Temperature: 310 K
Pressure: 1.0 bar
Protein force field: amber99sb-ildn
Water model: tip3p

Required before execution:
1. Manually validate receptor construct and binding pocket.
2. Convert the selected docking pose into a protein–DINP complex.
3. Generate DINP parameters with ACPYPE/GAFF (or an explicitly documented alternative).
4. Merge ligand topology with the GROMACS protein topology.
5. Energy minimization -> NVT -> NPT -> production MD.
6. Report RMSD, RMSF, radius of gyration, SASA, H-bonds, and binding free energy (MM-PBSA/GBSA).

Suggested commands (adapt paths/topologies after ligand parameterization):

    gmx pdb2gmx -f complex.pdb -o processed.gro -water tip3p -ff amber99sb-ildn
    gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt dodecahedron
    gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top
    gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
    gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
    gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr && gmx mdrun -deffnm em
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr && gmx mdrun -deffnm nvt
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr && gmx mdrun -deffnm npt
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr && gmx mdrun -deffnm md
