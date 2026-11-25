# USP26 protein: GROMACS setup, energy minimization and equilibration

> **Note:** this is the *standard* GROMACS workflow (solvate, neutralize, energy-minimize,
> NVT and NPT equilibration), following the public GROMACS lysozyme tutorial by
> Justin A. Lemkul: <http://www.mdtutorials.com/gmx/lysozyme/index.html>.
> It is not a novel method. The commands are listed here for transparency; the `.mdp`
> parameter files are **not** included in this repository (see README).
>
> Replace `USP26_model.pdb` with the actual USP26 structure file used (experimental
> structure or predicted model). The GROMACS version used is not recorded here.
>
> This document covers the steps up to and including NPT equilibration. A production
> MD run is not documented here.

```bash
source /usr/local/gromacs/bin/GMXRC

gmx --version

# Deleting H2O molecules (USP26_model.pdb = input USP26 structure)
grep -v HOH USP26_model.pdb > USP26_clean.pdb

# Preparing necessary files
# pdb2gmx asks interactively for the force field and water model; select the ones
# actually used. (CHARMM27 files are inspected below.)
gmx pdb2gmx -f USP26_clean.pdb -o USP26_processed.gro -ignh -ss

# Finding the force-field files if necessary
cd /usr/local/gromacs/share/gromacs/top
ls
cd charmm27.ff
ls
cd -

# Making a box
gmx editconf -f USP26_processed.gro -o box.gro -bt cubic -d 1.0 -c

# Adding solvent
gmx solvate -cp box.gro -cs spc216.gro -o sol.gro -p topol.top
tail topol.top

# Adding ions
gmx grompp -f ions.mdp -c sol.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o Ions-sol.gro -p topol.top -pname NA -nname CL -neutral
tail topol.top

# Energy minimization
# producing tpr file
gmx grompp -f minim.mdp -c Ions-sol.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
# checking system's energy (select "Potential" when prompted)
gmx energy -f em.edr -o Potential.xvg

# NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
nohup gmx mdrun -v -deffnm nvt &
tail nohup.out
# select "Temperature" when prompted
gmx energy -f nvt.edr -o Temperature.xvg

# NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# checking pressure (select "Pressure" when prompted)
gmx energy -f npt.edr -o Pressure.xvg
# checking density (select "Density" when prompted)
gmx energy -f npt.edr -o Density.xvg
```
