import sys,re
import numpy as np
from scipy import stats
import matplotlib as mpl 
mpl.use('agg')
import matplotlib.pyplot as plt

OAHOnly = False
OAMeOnly = False
CalcPairDiffs = False
CorrectOA = False
CorrectCB = False
WithUncert = True 
WithRep = True  

### Arguments: <Flags> <ExpType> <ExperimentFile> <CalculationFile1> [<CalculationFile2> ...]
### Note, files should have the same number of data points
calcfiles=[]
for arg in sys.argv:
  if arg == 'OAHOnly':
    OAHOnly = True
  if arg == 'OAMeOnly':
    OAMeOnly = True
  if arg == 'CalcPairDiffs':
    CalcPairDiffs = True
  if arg == 'CorrectOA':
    CorrectOA = True
  if arg == 'CorrectCB':
    CorrectCB = True
  if re.search(r'^(ka|dg|dh)$', arg):
    exptype = arg
  if re.search(r'\.txt$', arg):
    if re.search(r'Exp', arg):
      expfile = arg
    else:
      calcfiles.append(arg)

### Settings
R = 0.0019872036
T = 298.0
BootCyc = 100
MNames = ('Slope', 'Interc', 'R', 'R^2', 'RMSE', 'MSE', 'MUE', 'TAU')
Nm = len(MNames)

### Error Metric Calculations
def geterrormetrics(x, y):
  MTmp = np.zeros([Nm], np.float64)
  # Slope, Intercept, R
  MTmp[0], MTmp[1], MTmp[2], pval, stderr = stats.linregress(x,y)
  # R^2
  MTmp[3] = MTmp[2]**2
  # RMSE
  MTmp[4] = np.sqrt(np.mean(((y - x)**2)))
  # MSE
  MTmp[5] = np.mean((y - x))
  # MUE
  MTmp[6] = np.mean(np.absolute(y - x))
  # Tau
  MTmp[7], prob = stats.kendalltau(x,y)
  return (MTmp)


### BootStrapping Definition
def bootstrap(x, xsem, y, ysem):
  MBoot = np.zeros([Nm,BootCyc], np.float64)
  MVals = np.zeros([Nm], np.float64)
  MSEMs = np.zeros([Nm], np.float64)
  xtmp = np.zeros([len(x)], np.float64)
  ytmp = np.zeros([len(x)], np.float64)
  yfit = np.zeros([len(x)], np.float64)

  for b in range(BootCyc):
    for i in range(len(x)):

      # Sample with/without replacement?
      if WithRep:
        j = np.random.randint(len(x))
      else:
        j = i

      # Sampling Statistical Uncertainty
      if not WithUncert or xsem[j] == 0.0:
        xtmp[i] = x[j]
      else:
        xtmp[i] = np.random.normal(x[j],xsem[j])
      if not WithUncert or ysem[j] == 0.0:
        ytmp[i] = y[j]
      else:
        ytmp[i] = np.random.normal(y[j],ysem[j])

      MBoot[0:Nm,b] = geterrormetrics(xtmp, ytmp)
  
  for m in range(Nm):
    MVals[m]=np.mean(MBoot[m])
    MSEMs[m]=np.std(MBoot[m])

  return (MVals,MSEMs,MBoot)



### Load experimental file and place data in array
print expfile
with open(expfile, 'r') as expraw:
  explines = expraw.readlines()
exp = []
for line in explines:
  if not re.match(r'^\s*$', line):
    exp.append(line.rstrip().replace('\t', ' '))

### If "Only" flag, just do first six, or last six
if OAHOnly:
  exp = exp[0:6]
if OAMeOnly:
  exp = exp[6:12]
if OAHOnly and OAMeOnly:
  print "OAHOnly=True and OAMeOnly=True! Not compatible"
  exit()
N = len(exp)	# Number of data points
Np = (N-1)*N/2  # Number of data pairs

### Are these binding constants or free energy (or enthalpy)? Convert.
emean = np.zeros([N], np.float64)
esem = np.zeros([N], np.float64)
for i in range(len(exp)):
  cols = np.asarray(exp[i].split(), dtype=np.float64)
  if exptype == 'ka':       ### Assume binding constant. This could go wrong for very positive enthalpy
    dG = -R*T*np.log(cols)
    emean[i] = np.mean(dG)
    esem[i] = np.std(dG, ddof=1)/np.sqrt(len(dG))
  elif exptype == 'dg':
    emean[i] = cols[0]
    esem[i] = cols[1]
  elif exptype == 'dh': 
    emean[i] = np.mean(cols)/1000
    esem[i] = np.std(cols, ddof=1)/np.sqrt(len(cols))/1000
  else:
    print exptype, "... is not a valid experimental type"

### Calculate Experimental Pairwise Differences
if CalcPairDiffs:
  h = 0
  epmean = np.zeros([Np], np.float64)
  epsem = np.zeros([Np], np.float64)
  for i in range(len(exp)):
    for j in range(i+1, len(exp)):
      epmean[h] = emean[i] - emean[j]
      epsem[h] = np.sqrt( esem[i]**2 + esem[j]**2 )
      h += 1

### Read in Calculated data.  
### I should add a check to make sure number of data points is the same as experiment
Nc = len(calcfiles)
cmean = np.zeros([Nc,N], np.float64)
csem = np.zeros([Nc,N], np.float64)
if CalcPairDiffs:
  cdmean = np.zeros([Nc,Nd], np.float64)
  cdsem = np.zeros([Nc,Nd], np.float64)

RawMs = np.zeros([Nc,Nm], np.float64)
AllMBoot = np.zeros([Nc,Nm,BootCyc], np.float64)
AllMVals = np.zeros([Nc,Nm], np.float64)
AllMSEMs = np.zeros([Nc,Nm], np.float64)

### Print Metrics Column Headings
#print "%30s " % ("Submission"),
#for Name in MNames:
#  print "%6s " % (Name),
#print ""

nc = 0
for calcfile in calcfiles:
  calc = np.loadtxt(calcfile, np.float64)
  if OAHOnly:
    calc = calc[0:6]
  if OAMeOnly:
    calc = calc[6:12]
  for i in range(len(calc)):
    if np.isscalar(calc[i]) == True:      ### If scalar instead of array; ie, no SEM given.
      cmean[nc,i] = np.mean(calc[i])
      csem[nc,i] = 0.0
    else:                                 ### Assume Mean and SEM given
      cmean[nc,i] = calc[i,0]
      csem[nc,i] = calc[i,1]
    #print emean[i],esem[i],tmpcmean[i],tmpcsem[i]


  ### Correct data set with MSE
  if CorrectOA:
    cmean[nc,0:6] = cmean[nc,0:6] - (np.mean(cmean[nc,0:6]) - np.mean(emean[0:6]))
    if len(calc) == 12:
      cmean[nc,6:12] = cmean[nc,6:12]- (np.mean(cmean[nc,6:12]) - np.mean(emean[6:12]))
  if CorrectCB:
    cmean[nc,0:10] = cmean[nc,0:10] - (np.mean(cmean[nc,0:10]) - np.mean(emean[0:10]))

  if CalcPairDiffs:
    h = 0
    for i in range(len(calc)):
      for j in range(i+1, len(calc)):
        cpmean[nc,h] = cmean[nc,i] - cmean[nc,j]
        cpsem[nc,h] = np.sqrt( csem[nc,i]**2 + csem[nc,j]**2 )
        #print cpmean[nc,h]#,cdsem[nc,h]
        h += 1

  RawMs[nc] = geterrormetrics(emean, cmean[nc])
  AllMVals[nc],AllMSEMs[nc],AllMBoot[nc] = bootstrap(emean, esem, cmean[nc], csem[nc])
  #print "%30s " % (calcfile),
  #for m in range(Nm):
  #  print "%6.2f " % (RawMs[nc,m]),
  #print ""
  
  nc += 1



CalcNames=[]
for name in calcfiles:
  cols = name.split('.')
  CalcNames.append(cols[0])
  
Order = np.argsort(AllMVals[0:Nc,4])

Locs = np.arange(Nc)+1
plt.figure(1, figsize=(18,12), dpi=300)
plot2 = plt.boxplot([AllMBoot[nc,4] for nc in Order], widths=0.6, notch=True, sym='', positions=Locs)
ax = plt.axes()
ax.yaxis.grid(linestyle='--', which='major', color='lightgrey', alpha=0.7, zorder=0)
plt.setp(plot2['boxes'], color='b', linewidth=2.0) #,  facecolor='DarkMagenta')
plt.setp(plot2['whiskers'], color='b', linewidth=2.0, linestyle=':' )
plt.setp(plot2['caps'], color='b', linewidth=2.0)
plt.setp(plot2['medians'], color='None', linewidth=2.0)
plt.plot(Locs, [AllMVals[nc,4] for nc in Order], '_', color='OrangeRed', markersize=20, markeredgecolor='OrangeRed', markeredgewidth=3)
plt.plot(Locs, [RawMs[nc,4] for nc in Order], 'ko', markersize=8, markeredgecolor='k', markeredgewidth=1)
plt.xticks(Locs, [CalcNames[nc] for nc in Order], rotation=50, ha='right')
#plt.show()
plt.savefig('test.png', bbox_inches='tight')








#for nc in range(Nc):
#  print RawMs[nc]
#  print AllMVals[nc]

#fig = plt.figure(1, figsize=(9, 6))
#ax = fig.add_subplot(111)
#bp = ax.boxplot([AllMBoot[0:1,4], AllMBoot[1,4]])
#1fig.savefig('fig1.png', bbox_inches='tight')


