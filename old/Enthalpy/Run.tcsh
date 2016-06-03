#!/bin/tcsh
set echo
python ErrorAnalysis.py ExpEnthalpy10.txt APR_OPC-10.txt APR_TIP3P-10.txt Null.txt Dumb.txt Null2.txt > results.dat 
