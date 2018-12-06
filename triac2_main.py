from math import pi, pow
import matplotlib.pyplot as plt
import numpy as np

from triac2_utils import *
from InputData import *

# Dummy properties
inpf = 'inptriac2'
data = readdata(inpf)
diff = data['diff']
geod = data['geod']
print('data transport: ',data['diff'] )
print('------')
print('data geometry GEOD: ',data['geod'] )
print('------')
print('data geometry GEODI: ',data['geodi'] )
print('------')
print('data Kontaminasi HM: ',data['kont'] )
print('------')
print('data Recoil: ',data['recr'] )
print('------')
print('data Nukdat: ',data['nukdat'] )
print('------')
print('data Adsoc: ',data['adsoc'] )
print('------')
#
## part of ANFANG soubroutine (to be moved to separate subroutine/def after finished )
print('######')
print('ANFANG')
print('######')
#
ainv=data['nukdat']['ainv']
ukongp=data['kont']['ukongp']
aibe=ainv*ukongp
print('aibe : ',aibe)
dummyv=0.
dummyp=0.
dummyt=400.  #temp value
uez=50.      #temp value
npb=data['geodi']['npb']
print('npb: ', npb)
diff=data['diff']
di=[]
q=[]
difgp=[]
for i in range(npb):
    di.append(difpo(diff,dummyt,i))
    q.append(1.)
difgp.append(di[0])
dummy1=0.
dummy2=0.
#kud=kudif(npb,rb,dummy1, uez, di,q,zerfk, dummy2,cb,ges)

print('ukongp: ', ukongp)
#end of ANFANG
####
nd = 100
rd = np.zeros((nd+1), dtype=float)
td = np.zeros((nd+1), dtype=float)
tdix = np.zeros((nd+1), dtype=float)
rd[0] = 0
rd[1] = 5
td[0] = 100
for i in range (2, nd+1):
    rd[i] = rd[i-1] / 1.1

for i in range(0, nd):
    td[i+1] = td[i] * 0.98
td[nd] = 0
tdi = tuple(td)
for i in range(0, len(td)):
    tdix[i] = td[i]
dtid = 0.5
uezd = 50   # In graphite grain and in particle is set 50 CM^2/S
qd = np.ones((nd), dtype=float)
zerfkd = 10
t0d = 0
gesd = 0

# High --> Low (Diff. Const.)
did = np.ones((nd), dtype=float)
did[0] = 10

for i in range(0, nd - 1):
    did[i+1] = did[i] * 0.98 #  ((DID[I]**2)*0.08)

test= kudif(nd, rd, dtid, uezd, did, qd, zerfkd, t0d, td, gesd)
print('test : ', test)
td=test[0]
print('td : ', len(td), td)
# Plotting
l = np.arange(nd+1)

plt.plot(l, td, '--')
plt.title('Grafik Konsentrasi pada Tiap Layer (HI -> LO)')
plt.xlabel('Layer Ke - ')
plt.ylabel('Konsentrasi')
plt.legend(('Variasi 1'),loc = 0)
plt.grid(True)
plt.savefig('Concentration_HiLo.png',format ='png')
plt.show()

plt.show()

#input("\npress return to exit")

