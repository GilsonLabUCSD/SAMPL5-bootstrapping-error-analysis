#!/bin/tcsh
set echo
python ErrorAnalysis.py ExpSaltDep-OAH.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt CCSD\(T\)-neutral.txt DFT-charged.txt DFT-neutral.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt PERT-combo-OAH.txt Null.txt Dumb.txt Null2.txt > results.dat

