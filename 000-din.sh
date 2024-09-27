#!/bin/bash
proteina=PARP.pdb
ligando=GK4.itp
for arch in $proteina
do
  prot=$(echo $arch | sed 's/.pdb//g' )
  cp $prot.pdb old-p
  liga=$(echo $ligando | sed 's/.itp//g' )
  comp=$liga-$prot
  cp $comp.pdb old-c

gmx editconf -f $prot.pdb -o $prot.gro -resnr 1 -c -bt cubic -d 1
  gmx pdb2gmx -f $prot.gro -o $prot.gro -p $prot.top -i $prot.itp -ff oplsaa -his -water tip3p
 sed -i "/16008 16007 16003 16009/ s/$/ \n\n[ distance_restraints ]\n; ai   aj   type   index   type'      low     up1     up2     fac\n16014  319     1       0       1      0.0     0.15    0.28    1.0\n16014  373     1       1       1      0.0     0.15    0.28    1.0\n16014  833     1       2       1      0.0     0.15    0.23    1.0\n16014  880     1       3       1      0.0     0.15    0.28    1.0\n16015  4652    1       4       1      0.0     0.15    0.28    1.0\n16015  4692    1       5       1      0.0     0.15    0.28    1.0\n16015  4889    1       6       1      0.0     0.15    0.28    1.0\n16015  5032    1       7       1      0.0     0.15    0.28    1.0/g" $prot.top
  sed -i "/oplsaa.ff\/forcefield.itp/ s/$/ \n\#include \"$liga.itp\"/g" $prot.top
  sed -i "/Protein             1/ s/$/ \n$liga                 1/g" $prot.top
  gmx editconf -f $comp.pdb -o $comp.gro -c -bt cubic -d 1
  gmx solvate -cp $comp.gro -cs spc216.gro -o $comp.gro -p $prot.top
  gmx grompp -f min.mdp -c $comp.gro -p $prot.top -o $comp.tpr -po $comp.mdp -maxwarn 2
  echo SOL | gmx genion -s $comp.tpr -o $comp.gro -p $prot.top -neutral -conc 0.154004106
  gmx_d grompp -f min.mdp -c $comp.gro -p $prot.top -o $comp-min.tpr -po $comp-min.mdp
  gmx_d mdrun -nice 0 -v -s -pin on -pinoffset 0 -deffnm $comp-min
  
  ############################################

  gmx grompp -f nvt.mdp -c $comp-min.gro -p $prot.top -o $comp-nvt.tpr -po $comp-nvt.mdp
  gmx mdrun -gpu_id 0 -nice 0 -v -s -pin on -pinoffset 0 -nb gpu -deffnm $comp-nvt

  gmx grompp -f npt.mdp -c $comp-nvt.gro -p $prot.top -o $comp-npt.tpr -po $comp-npt.mdp
  gmx mdrun -gpu_id 0 -nice 0 -v -s -pin on -pinoffset 0 -nb gpu -deffnm $comp-npt
done
rm \#*
rm step*

######################################REPRENDER################################
#gmx convert-tpr -s ABC-PROT-npt.tpr -extend 50000 -o ABC-PROT-npt-2.tpr
#gmx mdrun -gpu_id 0 -nice 0 -v -deffnm ABC-PROT-npt -s ABC-PROT-npt-2.tpr -cpi ABC-PROT-npt.cpt -pin on -pinoffset 0 -nb gpu
######################################INDEX####################################
#gmx make_ndx -f ABC-PROT-npt.gro -o index.ndx
#####################################QUITAR PERIODICIDAD#######################
#gmx trjconv -f ABC-PROT-npt.xtc -s ABC-PROT-npt.tpr -ur compact -pbc res -center -n index.ndx -o ABC-NOPBC.xtc -dt 100
#####################################RMSD######################################
#gmx rms -s ABC-PROT-npt.tpr -f ABC-NOPBC.xtc -n index.ndx -o ABC-rmsd.xvg -tu ns
####################################RADIO DE GIRO##############################
#gmx gyrate -s ABC-PROT-npt.tpr -f ABC-NOPBC.xtc -n index.ndx -o ABC-gyrate.xvg
####################################RMSF#######################################
#gmx rmsf -f ABC-NOPBC.xtc -s ABC-PROT-npt.gro -o ABC-rmsf.xvg -res
###################################HBOND#######################################
#gmx hbond -f ABC-NOPBC.xtc -s ABC-PROT-npt.tpr -num ABC-hbond.xvg -tu ns
