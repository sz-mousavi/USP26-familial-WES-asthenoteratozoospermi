# USP26 Protein Minimization
# Commands used for the USP26 project (lab notebook)

source /usr/local/gromacs/bin/GMXRC
gmx --version
# Deleting H2O Molecules
grep -v HOH name.pdb > 1rex_clean.pdb

# ---------------------------------------------------------------------------
# Preparing-necessary-files
# ---------------------------------------------------------------------------
gmx pdb2gmx -f 1rex_clean.pdb -o 1rex.gro -ignh -ss

# ---------------------------------------------------------------------------
# finding ART file if necessary
# ---------------------------------------------------------------------------
cd source/user/local/gromacs/share/gromacs/top
ls
cd charmm27.ff
ls

# ---------------------------------------------------------------------------
# Making a box
# ---------------------------------------------------------------------------
gmx editconf -f 1rex.gro -o box.gro -bt cubic -d 1.0 -c

# ---------------------------------------------------------------------------
# Adding solvent
# ---------------------------------------------------------------------------
gmx solvate -cp box.gro -cs spc216.gro -o sol.gro -p topol.top
tail topol.top

# ---------------------------------------------------------------------------
# Adding ions
# ---------------------------------------------------------------------------
gmx grompp -f ions.mdp -c sol.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o Ions-sol.gro -p topol.top -pname NA -nname CL -neutral
tail topol.top

# ---------------------------------------------------------------------------
# Energy minimization
# ---------------------------------------------------------------------------
# producing tpr file
gmx grompp -f minim.mdp -c Ions-sol.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
# checking system's energy
gmx energy -f em.edr -o Potential.xvg

# ---------------------------------------------------------------------------
# NVT
# ---------------------------------------------------------------------------

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
nohup gmx mdrun -v -deffnm nvt &
tail nohup.out
gmx energy -f nvt.edr -o Temprature.xvg

# ---------------------------------------------------------------------------
# NPT
# ---------------------------------------------------------------------------
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt


# checking pressure
gmx energy -f npt.edr -o Pressure.xvg
# checking density
gmx energy -f npt.edr -o density.xvg
