#!/bin/bash

#Variables
protein = "PARP.pdb"
prot=$(basename $protein.pdb)
forcefield = "oplsaa_zn.ff"
solvate = true
ions = true

backupProt(){
  cp protein/${prot}.pdb $prot.pdb
  cp ${prot}.pdb ${prot}_B.pdb
}




