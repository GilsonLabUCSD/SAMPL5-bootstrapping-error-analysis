#!/bin/tcsh
set echo
python ErrorAnalysis.py ExpCBClip.txt BAR-ab-initio.txt BAR-dock.txt BAR-TI.txt BEDAM.txt MovTyp-1.txt MovTyp-2.txt SOMD-1.txt SOMD-2.txt SOMD-3.txt SOMD-4.txt TI-ab-initio.txt TI-dock.txt Null.txt Dumb.txt Null2.txt > results.dat
