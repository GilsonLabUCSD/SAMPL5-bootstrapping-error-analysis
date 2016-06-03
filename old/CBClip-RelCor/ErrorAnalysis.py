import numpy as np
from scipy import stats
import re
import sys 

### Arguments:  <ExperimentFile> <CalculationFile1> [<CalculationFile2> ...] 
### Note, files must have the same number of data points
expfile = sys.argv[1]
calcfiles = sys.argv[2:]

### Setup
R = 0.0019872036
T = 298.0
BootCyc = 100000

### Load Experimental Ka File and Prepare Arrays
print expfile
with open(expfile, 'r') as explines:
  expdata = explines.readlines()
exp = []
for line in expdata:
  if not re.match(r'^\s*$', line):
    exp.append(line.rstrip().replace('\t', ' '))
N = len(exp)
Nd = (len(exp)-1)*len(exp)/2
emean = np.zeros([N], np.float64)
esem = np.zeros([N], np.float64)
edmean = np.zeros([Nd], np.float64)
edsem = np.zeros([Nd], np.float64)

### Calculate Means and SEMs
for i in range(len(exp)):
  cols = np.asarray(exp[i].split(), dtype=np.float64)
  dG = -R*T*np.log(cols)
  emean[i] = np.mean(dG)
  #esem[i] = stats.sem(dG)  ### This is "sample standard deviation" divided by sqrt(n).  see http://stattrek.com/estimation/confidence-interval-mean.aspx?tutorial=stat#
  esem[i] = np.std(dG, ddof=1)/np.sqrt(len(dG))

### Alternate: Load Experimental dG +/- SEM file
#for i in range(len(exp)):
#  cols = np.asarray(exp[i].split(), dtype=np.float64)
#  #emean[i] = cols[0]
#  #esem[i] = cols[1]
#  emean[i] = np.mean(cols)/1000
#  #esem[i] = stats.sem(cols)/1000
#  esem[i] = np.std(cols, ddof=1)/np.sqrt(len(cols))/1000


### Calculate Relative Differences
#h = 0
#for i in range(len(exp)):
#  print emean[i],esem[i]
#  for j in range(i+1, len(exp)):
#    edmean[h] = emean[i] - emean[j]
#    edsem[h] = np.sqrt( esem[i]**2 + esem[j]**2 )
#    #print edmean[h]#,edsem[h]
#    h += 1

### Load Calculation Files and Run Statistics
Nc = len(calcfiles)
cmean = np.zeros([Nc,N], np.float64)
csem = np.zeros([Nc,N], np.float64)
cdmean = np.zeros([Nc,Nd], np.float64)
cdsem = np.zeros([Nc,Nd], np.float64)
etmp = np.zeros([N], np.float64)
ctmp = np.zeros([N], np.float64)
edtmp = np.zeros([Nd], np.float64)
cdtmp = np.zeros([Nd], np.float64)
R = np.zeros([BootCyc], np.float64)
R2 = np.zeros([BootCyc], np.float64)
RMSE = np.zeros([BootCyc], np.float64)
MSE = np.zeros([BootCyc], np.float64)
MUE = np.zeros([BootCyc], np.float64)
Slp = np.zeros([BootCyc], np.float64)
Int = np.zeros([BootCyc], np.float64)
Tau = np.zeros([BootCyc], np.float64)
dR = np.zeros([BootCyc], np.float64)
dR2 = np.zeros([BootCyc], np.float64)
dRMSE = np.zeros([BootCyc], np.float64)
dMSE = np.zeros([BootCyc], np.float64)
dMUE = np.zeros([BootCyc], np.float64)
dSlp = np.zeros([BootCyc], np.float64)
dInt = np.zeros([BootCyc], np.float64)
dTau = np.zeros([BootCyc], np.float64)

print  "     Name                 Slope       Interc.      R^2           RMSE          MSE          MUE          TAU"


nc = 0
for calcfile in calcfiles:
  calc = np.loadtxt(calcfile, np.float64)
  for i in range(len(calc)):
    if np.isscalar(calc[i]) == True:      ### No SEM given
      cmean[nc,i] = np.mean(calc[i])
      csem[nc,i] = 0.0
    else:                     ### Assume Mean and SEM given
      cmean[nc,i] = calc[i,0]
      csem[nc,i] = calc[i,1]
    print emean[i],esem[i],cmean[nc,i],csem[nc,i]

  ### Correct data set with MSE
  cmean[nc,0:10] = cmean[nc,0:10] - (np.mean(cmean[nc,0:10]) - np.mean(emean[0:10]))
#  if len(calc) == 12:
#    cmean[nc,6:12] = cmean[nc,6:12]- (np.mean(cmean[nc,6:12]) - np.mean(emean[6:12]))

  for i in range(len(calc)):
    print cmean[nc,i],csem[nc,i]


#  h = 0
#  for i in range(len(calc)):
#    for j in range(i+1, len(calc)):
#      cdmean[nc,h] = cmean[nc,i] - cmean[nc,j]
#      cdsem[nc,h] = np.sqrt( csem[nc,i]**2 + csem[nc,j]**2 )
#      #print cdmean[nc,h]#,cdsem[nc,h]
#      h += 1

  for b in range(BootCyc):
    for i in range(N):
      if esem[i] == 0.0:
        etmp[i] = emean[i]
      else:
        etmp[i] = np.random.normal(emean[i],esem[i])
      if csem[nc,i] == 0.0:
        ctmp[i] = cmean[nc,i]
      else:
        ctmp[i] = np.random.normal(cmean[nc,i],csem[nc,i])
    Slp[b], Int[b], R[b], pval, stderr = stats.linregress(etmp,ctmp)
    R2[b] = R[b]**2
    RMSE[b] = np.sqrt(np.mean(((ctmp - etmp) ** 2)))
    MSE[b] = np.mean((ctmp - etmp))
    MUE[b] = np.mean(np.absolute(ctmp - etmp))
    Tau[b], prob = stats.kendalltau(etmp,ctmp)

#    for i in range(Nd):
#      if edsem[i] == 0.0:
#        edtmp[i] = edmean[i]
#      else:
#        edtmp[i] = np.random.normal(edmean[i],edsem[i])
#      if cdsem[nc,i] == 0.0:
#        cdtmp[i] = cdmean[nc,i]
#      else:
#        cdtmp[i] = np.random.normal(cdmean[nc,i],cdsem[nc,i])
#    dSlp[b], dInt[b], dR[b], pval, stderr = stats.linregress(edtmp,cdtmp)
#    dR2[b] = dR[b]**2
#    dRMSE[b] = np.sqrt(np.mean(((cdtmp - edtmp) ** 2)))
#    dMSE[b] = np.mean((cdtmp - edtmp))
#    dMUE[b] = np.mean(np.absolute(cdtmp - edtmp))
#    dTau[b], prob = stats.kendalltau(edtmp,cdtmp)
  print "%30s  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f" % ( calcfile, np.mean(Slp), np.std(Slp), np.mean(Int), np.std(Int), np.mean(R2), np.std(R2), np.mean(RMSE), np.std(RMSE), np.mean(MSE), np.std(MSE), np.mean(MUE), np.std(MUE), np.mean(Tau), np.std(Tau) )
  #print "%30s  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f  %6.2f %6.2f" % ( calcfile, np.mean(dSlp), np.std(dSlp), np.mean(dInt), np.std(dInt), np.mean(dR2), np.std(dR2), np.mean(dRMSE), np.std(dRMSE), np.mean(dMSE), np.std(dMSE), np.mean(dMUE), np.std(dMUE), np.mean(dTau), np.std(dTau) )
  nc += 1

 
