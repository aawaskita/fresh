import numpy as np

geod={}
geod['rbe']=[]
geod['rcp']=[]
geod['rkorn']=0
geod['rb']=[]
geod['rp']=np.zeros(200).tolist() """defined as 200 element array in padiff sub routine"""
geod['rko']=[]
geod['pzahl0']=0

geodi={}
geodi['nrbe']=0
geodi['nrp']=0
geodi['npb']=0
geodi['npp']=0
geodi['nko']=0
geodi['nrb']=np.zeros(5).tolist() """defined as 5 element array in padiff, determined by nrbe otherwise"""
geodi['nrc']=np.zeros(5).tolist() """defined as 5 element array in padiff, determined by nrp otherwise"""

adsoc={}
adsoc['aci']=0
adsoc['bci']=0
adsoc['cgrenc']=0
adsoc['tgrenc']=0
adsoc['anci']=0
adsoc['bnci']=0
adsoc['enci']=0
adsoc['fnci']=0
adsoc['ifadc']=0

cpbru={}
cpbru['pbra']=0
cpbru['pzahli']=0
cpbru['pzahld']=np.zeros(10).tolist()
cpbru['freik']=0
cpbru['frcp']=0
cpbru['frb']=0
cpbru['frgk']=0

diff={}
diff['dog']=0
diff['akg']=0
diff['dok']=0
diff['akk']=0
diff['dop']=[]
diff['akp']=[]
diff['fog']=0
diff['fogk']=0
diff['fopk']=0
diff['fopsic']=0
diff['fodpk']=0

frerat={}
frerat['freii']=0
frerat['freid']=0
frerat['frecpg']=0
frerat['freiko']=0
frerat['freibe']=0

invent={}
invent['aicp']=[] """jumlah elemennya ditentukan nrp"""
invent['aicpi']=0
invent['aicpd']=np.zeros(10).tolist()
invent['aicpk']=0
invent['aicpdg']=[]
invent['aicpg']=0
invent['aigk']=0

inven1={}
inven1['aigp']=0
inven1['aifr']=0

kont={}
kont['ukongk']=0
kont['ukongp']=0
kont['ukontp']=[]

konz={}
konz['cb']=np.zeros(200).tolist()
konz['ck']=np.zeros(200).tolist()
konz['cp']=np.zeros(200).tolist()
konz['cpb']=np.zeros(shape=(50,10),dtype=float)

nukdat={}
nukdat['ainv']=0
nukdat['zerfk']=0

outst={}
outst['iopt']=np.zeros(20).tolist()
outst['ioptpl']=0
outst['ifplot']=0
outst['ifbiso']=0
outst['ifsick']=0
outst['ngnr']=0
outst['ifplt1']=0
outst['ifplt2']=0

parrel={}
parrel['reld']=[]
parrel['reldi']=0
parrel['relii']=0
parrel['ibr']=0

pp={}
pp['x1']=[]
pp['y1']=[]
pp['y2']=[]
pp['y3']=[]
pp['y4']=[]
pp['y5']=[]
pp['ilo']=0

quell={}
quell['zeit0']=0
quell['qgesam']=0.0
quell['qakt']=0
quell['q1']=0.0

real1={}
real1['zeit0d']=0
real1['xn0']=0
real1['zeitpr']=0
real1['gamma']=0

stopp={}
stopp['ifstop']=0

randb={}
ranb['stoffu']=0
ranb['alpha']=0
ranb['c00k']=0

recfr={}
recfr['rcfrpk']=0
recfr['rcfrp']=0
recfr['rcfrdf']=0
recfr['rcfra']=0
recfr['rcfrbe']=0

recr={}
recr['recker']=0
recr['recpyc']=0
recr['recgra']=0

relea={}
relea['relint']=0
relea['reldef']=0
relea['relk']=0
relea['relcpg']=0
relea['relgk']=0
relea['relbe']=0

relrec={}
relrec['relrpk']=0
relrec['relrp']=0
relrec['relrdf']=0
relrec['relra']=0
relrec['relrbe']=0

tradat={}
tradat['temper']=0
tradat['difgp']=np.zeros(5).tolist()
tradat['difgk']=0
tradat['difkcp']=np.zeros(5).tolist()

ueber=''
title=''

zeits={}
zeits['te']=[]
zeits['dt']=[]
zeits['dtout']=[]

zeit2={}
zeit2['zeit']=0
zeit2['dzeit']=0

rio={}
rio['nout']=0
rio['nin']=0
rio['nfrk']=0
rio['nfrkm']=0

sic={}
sic['alsic']=0.0
sic['ds0']=0
sic['dkor']=0
sic['dsic']=0
sic['dif']=0
sic['isic']=0
sic['msic']=0
sic['msiges']=0

adsoc={}
adsoc['aci']=0
adsoc['bci']=0
adsoc['cgrenc']=0
adsoc['tgrenc']=0
adsoc['anci']=0
adsoc['bnci']=0
adsoc['enci']=0
adsoc['fnci']=0
adsoc['ifadc']=0

int1={}
int1['lrcp']=0
int1['lncp']=0
int1['lzs']=0

int2={}
int2['ircp']=0
int2['izs']=0
int2['ns']=0

"""Variables that defined as DIMENSIONS in fortran. Most of them were returned back to list after declaring using numpy, while some of them remain in numpy format"""

"""a=np.zeros(200).tolist() declared as different number of elements"""
abw=np.zeros(60).tolist()
acp=np.zeros(60).tolist()
afrei=np.zeros(60).tolist()
agk=np.zeros(60).tolist()
agp=np.zeros(60).tolist()
akern=np.zeros(60).tolist()
"""b=np.zeros(200).tolist() declared as different number of elements"""
betas=np.zeros(200).tolist()
bil=np.zeros(60).tolist()
"""c=np.zeros(200).tolist() in kudif subroutine, otherwise it was declared as 5 element array"""
cbruch=np.zeros(50).tolist()
cpdef=np.zeros(60).tolist()
d=np.zeros(200).tolist()
di=np.zeros(200).tolist()
did=np.zeros(200).tolist()
difdcp=np.zeros(5).tolist()
difg=np.zeros(shape=(5,60), dtype=float)
difk=np.zeros(60).tolist()
difp=np.zeros(shape=(5,60), dtype=float)
dik=np.zeros(200).tolist()
dip=np.zeros(200).tolist()
frbe=np.zeros(60).tolist()
freiby=np.zeros(10).tolist()
freidy=np.zeros(10).tolist()
frk=np.zeros(60).tolist()
frp=np.zeros(60).tolist()
isp=np.zeros(60).tolist()
iopt1=np.zeros(20).tolist()
q=np.zeros(200).tolist()
q11=np.zeros(60).tolist()
qa=np.zeros(60).tolist()
qk=np.zeros(200).tolist()
qp=np.zeros(200).tolist()
r=np.zeros(200).tolist()
"""sp=np.zeros(shape=(5,60),dtype=float) in out18 sub routine, it was declared as 3,60 2D element matrix"""
t=np.zeros(200).tolist()
tem=np.zeros(60).tolist()
v=np.zeros(200).tolist()
zeit1=np.zeros(60).tolist()

