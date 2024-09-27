#!/bin/bash

# Variables
proteina="PARP-1.pdb"  # Nombre del archivo .pdb de la proteína
ligando="GK4.itp"       # Nombre del archivo .itp del ligando
toppar="Metalloprotein/toppar"  # Directorio donde se encuentran los archivos .itp
forcefield="oplsaa"     # Campo de fuerza
solvate=true            # Para saber si se debe solvatizar
ions=true          # Para saber si se deben añadir iones

# Crear un nombre base para la proteína
prot=$(basename "$proteina" .pdb)

# Copiar la proteína a una copia temporal
cp "Metalloprotein/protein/$proteina" "$prot.pdb"

# Preparar la estructura de la proteína
gmx editconf -f "$prot.pdb" -o "$prot.gro" -resnr 1 -c -bt cubic -d 1
gmx pdb2gmx -f "$prot.gro" -o "$prot.gro" -p "$prot.top" -i "$toppar/$prot.itp" -ff "$forcefield" -his -water tip3p

# Modificar el archivo .top para incluir el ligando
liga=$(basename "$ligando" .itp)
sed -i "/oplsaa.ff\/forcefield.itp/ s/$/ \n#include \"$toppar/$ligando\"/g" "$prot.top"
sed -i "/Protein             1/ s/$/ \n$liga                 1/g" "$prot.top"

# Solvatar el sistema si es necesario
if [ "$solvate" = true ]; then
  gmx solvate -cp "$prot.gro" -cs spc216.gro -o "$prot-solv.gro" -p "$prot.top"
fi

# Preparación para la minimización de energía
gmx grompp -f "Metalloprotein/mdp/min.mdp" -c "$prot-solv.gro" -p "$prot.top" -o "$prot.tpr" -po "$prot.mdp" -maxwarn 2

# Añadir iones para neutralizar el sistema si se ha solvado
if [ "$solvate" = true ] && [ "$ions" = true ]; then
  echo "Añadiendo iones..."
  echo SOL | gmx genion -s "$prot.tpr" -o "$prot-solv.gro" -p "$prot.top" -neutral -conc 0.154004106
  gmx grompp -f "Metalloprotein/mdp/min.mdp" -c "$prot-solv.gro" -p "$prot.top" -o "$prot-min.tpr" -po "$prot-min.mdp" -maxwarn 2
fi

# Ejecutar la minimización de energía
gmx mdrun -v -s "$prot-min.tpr" -deffnm "$prot-min"

############################################
# Simulación NVT
gmx grompp -f "Metalloprotein/mdp/nvt.mdp" -c "$prot-min.gro" -p "$prot.top" -o "$prot-nvt.tpr" -po "$prot-nvt.mdp"
gmx mdrun -gpu_id 0 -nice 0 -v -s "$prot-nvt.tpr" -deffnm "$prot-nvt"

# Simulación NPT
gmx grompp -f "Metalloprotein/mdp/npt.mdp" -c "$prot-nvt.gro" -p "$prot.top" -o "$prot-npt.tpr" -po "$prot-npt.mdp"
gmx mdrun -gpu_id 0 -nice 0 -v -s "$prot-npt.tpr" -deffnm "$prot-npt"

# Limpieza de archivos temporales
rm \#*
rm step*
