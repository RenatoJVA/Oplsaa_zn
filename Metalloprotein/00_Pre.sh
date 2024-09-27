#!/bin/bash

# Variables
protein="PARP-1.pdb"  # Sin espacios alrededor del igual
prot=$(basename "$protein" .pdb)  # Obtener el nombre de la proteína sin la extensión
forcefield="oplsaa_zn.ff"  # Campo de fuerza
solvate=true  # Opción para solvatar el sistema
ions=true     # Opción para añadir iones
toppar="Metalloprotein/toppar"  # Directorio para los archivos .itp
input_dir="Metalloprotein/protein"  # Directorio donde se encuentra el archivo .pdb
mdp_dir="Metalloprotein/mdp"  # Directorio donde se encuentran los archivos .mdp

# Crear el directorio toppar si no existe
mkdir -p "$toppar"

# Verificar si el archivo .pdb existe
if [ ! -f "$input_dir/$protein" ]; then
  echo "Error: No se encontró el archivo $input_dir/$protein"
  exit 1
fi

# Copiar el archivo PDB al directorio actual y crear una copia
cp "$input_dir/$protein" .  # Copia el archivo PDB a la ubicación actual
cp "$protein" "${prot}_B.pdb"  # Crear una copia del archivo

# Crear archivo .gro a partir del archivo .pdb
gmx editconf -f "$protein" -o "$prot.gro" -resnr 1 -c -bt cubic -d 1

# Generar topología, incluir campo de fuerza local y guardar .itp en toppar
gmx pdb2gmx -f "$prot.gro" -o "$prot.gro" -p "$prot.top" -i "$toppar/$prot.itp" -ff "$forcefield" -his -water tip3p

# Incluir archivos .itp adicionales en el archivo .top
{
  echo "#include \"$toppar/tip3p.itp\""
  echo "#include \"$toppar/ions.itp\""
  echo "#include \"$toppar/zinc.itp\""  # Incluir Zn
  echo "#include \"$toppar/oxt.itp\""    # Incluir OXT
} >> "$prot.top"

# Solvatar el sistema si es necesario
if [ "$solvate" = true ]; then
  gmx solvate -cp "$prot.gro" -cs spc216.gro -o "$prot.gro" -p "$prot.top"
fi

# Preparación para la minimización de energía antes de añadir iones
gmx grompp -f "$mdp_dir/min.mdp" -c "$prot.gro" -p "$prot.top" -o "$prot.tpr" -po "$prot.mdp" -maxwarn 2

# Neutralizar el sistema con iones si se solvató
if [ "$solvate" = true ] && [ "$ions" = true ]; then
  echo "Añadiendo iones..."
  echo SOL | gmx genion -s "$prot.tpr" -o "$prot.gro" -p "$prot.top" -neutral -conc 0.154004106
  gmx grompp -f "$mdp_dir/min.mdp" -c "$prot.gro" -p "$prot.top" -o "$prot-min.tpr" -po "$prot-min.mdp" -maxwarn 2
fi

# Ejecutar la minimización de energía
gmx mdrun -v -s "$prot-min.tpr" -deffnm "$prot-min"

# Copiar los archivos .itp relevantes al directorio toppar (esto es opcional si ya están ahí)
cp "Metalloprotein/oplsaa_zn.ff/"*.itp "$toppar/"

############################################
# Simulación NVT
gmx grompp -f "$mdp_dir/nvt.mdp" -c "$prot-min.gro" -p "$prot.top" -o "$prot-nvt.tpr" -po "$prot-nvt.mdp"
gmx mdrun -gpu_id 0 -nice 0 -v -s -pin on -pinoffset 0 -nb gpu -deffnm "$prot-nvt"

# Simulación NPT
gmx grompp -f "$mdp_dir/npt.mdp" -c "$prot-nvt.gro" -p "$prot.top" -o "$prot-npt.tpr" -po "$prot-npt.mdp"
gmx mdrun -gpu_id 0 -nice 0 -v -s -pin on -pinoffset 0 -nb gpu -deffnm "$prot-npt"

# Limpieza de archivos temporales
rm \#*
rm step*
