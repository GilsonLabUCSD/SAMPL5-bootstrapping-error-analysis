#!/bin/tcsh

set echo

python ../blah.py ka ../Exp/OAAllAvg.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt Null.txt Null2.txt | egrep 'txt|Sub' | sed 's/-nan/nan/g' > & AllAvg-Abs.dat

python ../blah.py ka ../Exp/OASaltDep.txt APR-OPC.txt APR-TIP3P.txt BEDAM.txt Metadynamics.txt MMPBSA-GAFF.txt MovTyp-1.txt MovTyp-2.txt PERT-bound-c.txt PERT-bound.txt PERT-hrex-c1.txt PERT-hrex-c2.txt PERT-hrex-c.txt PERT-hrex.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt Null.txt Null2.txt | egrep 'txt|Sub' | sed 's/-nan/nan/g' > & SaltDep-Abs.dat



