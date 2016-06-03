#!/bin/tcsh
set echo
python ErrorAnalysis.py ExpSaltDep.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt > results.dat
python ErrorAnalysis.py ExpSaltDep-OAH.txt PERT-combo-OAH.txt >> results.dat
python ErrorAnalysis.py ExpSaltDep.txt Null.txt >> results.dat
python ErrorAnalysis.py ExpSaltDep.txt Dumb.txt >> results.dat
python ErrorAnalysis.py ExpSaltDep.txt Null2.txt >> results.dat
