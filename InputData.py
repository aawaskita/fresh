import sys, math, re
import globalVar

def readdata(namafile):
	f=open(namafile,"r")
	"""Start Reading Input File"""
	
	a=f.readline()
	b=a.split('\t')
	c=b[0].split(' ')
	globalVar.outst['ngnr']=int(c[3])
	
	title=''
	for i in range(4):
		title=title+f.readline()
	globalVar.title=title
		
	te=[]
	temp1=f.readline()
	temp2=temp1.split()
	item=len(temp2)
	for i in range(item):
		te.append(float(temp2[i]))
	globalVar.zeits['te']=te
	
	dt=[]
	temp1=f.readline()
	temp2=temp1.split()
	for i in range(item):
		dt.append(float(temp2[i]))
	globalVar.zeits['dt']=dt
	
	dtout=[]
	temp1=f.readline()
	temp2=temp1.split()
	item=len(temp2)
	for i in range(item):
		dtout.append(float(temp2[i]))
	globalVar.zeits['dtout']=dtout
    
	uk=[]
	temp1=f.readline()
	temp2=temp1.split()
	temp3=int(temp2[0])
	uk.append(temp3)
	temp3=int(temp2[1])
	uk.append(temp3)
	temp3=int(temp2[2])
	uk.append(temp3)
	temp3=int(temp2[3])
	uk.append(temp3)
	temp3=int(temp2[4])
	uk.append(temp3)
	
	temp1=f.readline()
	temp2=temp1.split()
	temp3=int(temp2[0])
	uk.append(temp3)
	temp3=int(temp2[1])
	uk.append(temp3)
	
	temp1=f.readline()
	temp2=temp1.split()
	temp3=int(temp2[0])
	uk.append(temp3)
	temp3=int(temp2[1])
	uk.append(temp3)
	temp3=int(temp2[2])
	uk.append(temp3)
	temp3=int(temp2[3])
	uk.append(temp3)
	
	rbe=[]                      
	temp1=f.readline()
	temp2=temp1.split()
	item=len(temp2)
	for i in range(item):
		rbe.append(float(temp2[i]))
	"""outer radius of i-th zone in the in the fuel element (or pebble fuel)."""
	globalVar.geod['rbe']=rbe
	
	nrb=[]
	temp1=f.readline()
	temp2=temp1.split()
	item=len(temp2)
	npb=0
	for i in range(item):
		nrb.append(int(temp2[i]))
		npb=npb+int(temp2[i])
	"""number of mesh in each zone"""
	globalVar.geodi['nrb']=nrb
	
	"""number of total zones in the fuel element"""
	globalVar.geodi['nrbe']=item
	
	"""number of total meshes in fuel element (adding meshes from all zones)"""
	globalVar.geodi['npb']=npb
    
	temp1=f.readline()
	temp2=temp1.split()
	
	"""number of particles in fuel sphere"""
	globalVar.geod['pzahl0']=float(temp2[0])
	
	"""normal operation/irradiation time"""
	globalVar.real1['zeit0d']=float(temp2[1])
	
	"""total inventory of the fission product species considered"""
	globalVar.real1['xn0']=0
	
	"""time after begin of the accident/heating phase."""
	globalVar.real1['zeitpr']=float(temp2[2])
	
	"""fast neutron fluence at begin of the accident/heating phase"""
	globalVar.real1['gamma']=float(temp2[3])
    
	rcp=[]            
	"""first rcp, rcp[0] is the outer diam. of fuel kernel."""              
	temp1=f.readline()
	temp2=temp1.split()
	item=len(temp2)
	for i in range(item):
		rcp.append(float(temp2[i]))
		
	"""outer radius of i-th zone in the in the coated particle."""
	globalVar.geod['rcp']=rcp
	
	temp1=f.readline()
	temp2=temp1.split()
	"""radius of grahite grain (graphite pore is the grain boundary)"""
	globalVar.geod['rkorn']=float(temp2[0])
	
	nrc=[]
	temp1=f.readline()
	temp2=temp1.split()
	item=len(temp2)
	npp=0
	for i in range(item):
		nrc.append(int(temp2[i]))
		npp=npp+int(temp2[i])
	"""number of mesh in each zone"""
	globalVar.geodi['nrc']=nrc
	
	"""number of zones in fuel particle"""
	globalVar.geodi['nrp']=item
	
	"""number of total mesh in fuel particle (adding meshes in all zones)"""
	globalVar.geodi['npp']=npp
	
	temp1=f.readline()
	temp2=temp1.split()
	globalVar.geodi['nk0']=int(temp2[0])

	"""Geometry data calculated from raw input zone radius in fuel element : array 'rb'"""
	rb=[]
	"""it's the center of the fuel element, as first rb array."""
	rb.append(0.0)
	
	"""it was previously set to loop in a range from 1 to npb+1, but python start the index from 0
	warning: current algorithm is only suitable for nrbe=2 -> should we apply if to check the value of nrbe?"""
	for i in range(npb):        
		if i<=nrb[0]:
			mesh=rbe[0]/nrb[0]
			rb.append(rb[i-1]+mesh)
		else:
			mesh=(rbe[1]-rbe[0])/nrb[1]
			"""print(i, rb[i-1], mesh)"""
			rb.append(rb[i-1]+mesh)

	#print('-----')
	#print('radius of zone in fuel element: ', len(rb),rb)
	globalVar.geod['rb']=rb
	#
	# zone radius in fuel particle : array 'rp'
	mesh=[]
	accitem=0
	acc=[]
	nrp=globalVar.geodi['nrp']
	"""for i in range(nrp):"""
	for i in range(5):
		if i==0:
			mesh.append(rcp[0]/nrc[0])
			accitem=accitem+nrc[0]
			acc.append(accitem)
		else:
			mesh.append((rcp[i]-rcp[i-1])/nrc[i])
			accitem=accitem+nrc[0]
			acc.append(accitem)

	"""it's the center of the fuel particle, as first rp array
	warning: current algorithm is only suitable for nrbe=2 -> should we apply if to check the value of nrbe?"""
	rp=[]
	rp.append(0.0)                  
	temp=0.0
	for i in range(nrp):        
		for j in range(nrc[i]):
			rp.append(temp+mesh[i])
			"""temp value after the inner looping remain the same and added to the next level of the outer loop"""
			temp=temp+mesh[i]
	globalVar.geod['rp']=rp
	
	"""zone radius in graphite grain"""
	rk0=[]
	rk0.append(0.0)
	nk0=globalVar.geodi['nk0']
	rkorn=globalVar.geod['rkorn']
	for i in range(1,nk0+1):
		rk0.append(rk0[i-1]+(rkorn/nk0))

	globalVar.geod['rk0']=rk0
       
	"""Heavy Metal Contamination Data"""
	ukontp=[]
	temp1=f.readline()
	temp2=temp1.split()
	globalVar.kont['ukongk']=float(temp2[0])
	globalVar.kont['ukongp']=float(temp2[1])
	unkernel=1.-(float(temp2[0])+float(temp2[1])+float(temp2[2])+float(temp2[3])+float(temp2[4])+float(temp2[5]))
	ukontp.append(unkernel)
	ukontp.append(float(temp2[2]))
	ukontp.append(float(temp2[3]))
	ukontp.append(float(temp2[4]))
	ukontp.append(float(temp2[5]))
	globalVar.kont['ukontp']=ukontp

	"""Transport Data"""
	temp1=f.readline()
	temp2=temp1.split()
	temp3=len(temp2)
	globalVar.diff['d0g']=float(temp2[0])
	globalVar.diff['akg']=float(temp2[1])
	globalVar.diff['f0g']=float(temp2[2])
	globalVar.diff['f0pk']=float(temp2[3])
	globalVar.diff['f0psic']=float(temp2[4])
	globalVar.diff['f0pdpk']=float(temp2[5])
		
	temp1=f.readline()
	temp2=temp1.split()
	globalVar.diff['d0k']=float(temp2[0])
	globalVar.diff['akk']=float(temp2[1])
	globalVar.diff['f0kg']=float(temp2[2])
    
    	
	d0p=[]
	temp1=f.readline()
	temp2=temp1.split()
	temp3=float(temp2[0])
	d0p.append(temp3)
	temp3=float(temp2[1])
	d0p.append(temp3)
	temp3=float(temp2[2])
	d0p.append(temp3)
	temp3=float(temp2[3])
	d0p.append(temp3)
	temp3=float(temp2[4])
	d0p.append(temp3)
	globalVar.diff['d0p']=d0p
	
	akp=[]
	temp1=f.readline()
	temp2=temp1.split()
	temp3=float(temp2[0])
	akp.append(temp3)
	temp3=float(temp2[1])
	akp.append(temp3)
	temp3=float(temp2[2])
	akp.append(temp3)
	temp3=float(temp2[3])
	akp.append(temp3)
	temp3=float(temp2[4])
	akp.append(temp3)
	globalVar.diff['akp']=akp
	
	temp1=f.readline()
	temp2=temp1.split()
	globalVar.recr['recker']=float(temp2[0])
	globalVar.recr['recpyc']=float(temp2[1])
	globalVar.recr['recgra']=float(temp2[2])
    
	temp1=f.readline()
	temp2=temp1.split()
	globalVar.adsoc['aci']=float(temp2[0])
	globalVar.adsoc['bci']=float(temp2[1])
	globalVar.adsoc['cgrenc']=float(temp2[2])
	globalVar.adsoc['ifadc']=float(temp2[3])
	globalVar.adsoc['tgrenc']=0.0
	
	temp1=f.readline()
	temp2=temp1.split()
	globalVar.adsoc['anci']=float(temp2[0])
	globalVar.adsoc['bnci']=float(temp2[1])
	globalVar.adsoc['enci']=float(temp2[2])
	globalVar.adsoc['fnci']=float(temp2[3])

	temp1=f.readline()
	temp2=temp1.split()
	globalVar.nukdat['ainv']=float(temp2[0])
	if globalVar.outst['ngnr']==1:
		globalVar.nukdat['zerfk']=7.29E-10
	elif globalVar.outst['ngnr']==2:
		globalVar.nukdat['zerfk']=7.712E-10
	elif globalVar.outst['ngnr']==3:
		globalVar.nukdat['zerfk']=3.21E-8
	elif globalVar.outst['ngnr']==4:
		globalVar.nukdat['zerfk']=1.0E-6

	temp1=f.readline()      #to read the 'Irradiation Temp. History line
	temp1=f.readline()
	temp2=temp1.split()
	ni0=int(temp2[0])               # ni0 : number of pair of data in normal operation
	globalVar.itemp['ni0']=ni0
    
	for i in range(ni0):
		temp1=f.readline()
		temp2=temp1.split()
		globalVar.ttemp['azeit0'].append(float(temp2[1]))
		globalVar.ttemp['atemp0'].append(float(temp2[2]))
	"""globalVar.ttemp['azeit0']=azeit0          # units in [days]
	globalVar.ttemp['atemp0']=atemp0          # units in degC"""

	temp1=f.readline()      #to read the 'Accident Temp. History line
	temp1=f.readline()
	temp2=temp1.split()
	ni=int(temp2[0])
	globalVar.itemp['ni']=ni
    
	azeit=[]
	atemp=[]
	for i in range(ni):
		temp1=f.readline()
		temp2=temp1.split()
		azeit.append(float(temp2[1]))
		atemp.append(float(temp2[2]))
	globalVar.ttemp['azeit']=azeit                #units in [hours]
	globalVar.ttemp['atemp']=atemp                #units in degC

	alsic=0                 #corroded length of SiC layer in cm
	ds0=rcp[3]              #thickness of SiC layer , taken from rcp input array.
	dkor=0.0                # ?
	isic=nrc[3]             # number of SiC-layer
	dsic=ds0/isic           # thickness of a spherical shell in the SiC-layer.
	dif=0.0
	msic=0                  #number of corroded sphere shells in SiC
	msiges=0                # 
	globalVar.sic['alsic']=alsic
	globalVar.sic['ds0']=ds0
	globalVar.sic['dkor']=dkor
	globalVar.sic['dsic']=dsic
	globalVar.sic['dif']=dif
	globalVar.sic['isic']=isic
	globalVar.sic['msic']=msic
	globalVar.sic['msiges']=msiges

	"""fp-inventory in single layers"""
	aicp=[]
	for i in range(5):
		aicp.append(0.0)
	globalVar.invent['aicp']=aicp

	"""fp-inventory in intact coated particle (CP)"""
	aicpi=0.0 
	globalVar.invent['aicpi']=aicpi

	"""sp-inventory in one defective CP (kernel)"""
	aicpd=[]
	for i in range(5):
		aicpd.append(0.0)
	globalVar.invent['aicpd']=aicpd

	"""fp-inventory in kernels of all coated particle (CP)"""
	aicpk=0.0
	globalVar.invent['aicpk']=aicpk
    
	"""fp-inventory in kernels of all defective coated particle (CP)"""
	aicpdg=0.0
	globalVar.invent['aicpdg']=aicpdg

	"""fp-inventory in all coated particles (CPs)"""
	aicpg=0.0
	globalVar.invent['aicpg']=aicpg

	"""fp-inventory in dern graphitekoernern"""
	aigk=0.0
	globalVar.invent['aigk']=aigk

	"""fp-inventory in graphite-pores"""
	aigp=0.0
	globalVar.inven1['aigp']=aigp

	"""fp-amount outside the fuel element"""
	aifr=0.0
	globalVar.inven1['aifr']=aifr 

	"""fp-concentration in graphite pores"""
	cb=[]
	icb=200
	for i in range(icb):
		cb.append(0.0)
	globalVar.konz['cb']=cb 

	"""fp-concentration in graphite grains"""
	ck=[]
	ick=200
	for i in range(ick):
		ck.append(0.0)
	globalVar.konz['ck']=ck

	"""fp-concentration in graphite particles"""
	cp=[]
	icp=200
	for i in range(icp):
		cp.append(0.0)
	globalVar.konz['cp']=cp

	"""fp-concentration in graphite grains"""
	cpb=[]
	icpb=200
	for i in range(icpb):
		cpb.append(0.0)
	globalVar.konz['cpb']=cpb
 
	"""temperature"""
	temper=0.0
	globalVar.tradat['temper']=temper 
    
	"""fp-concentration in graphite grains"""
	difgp=[]
	idifgp=5
	for i in range(idifgp):
		difgp.append(0.0)
	globalVar.tradat['difgp']=difgp

	"""fp-amount outside the fuel element"""
	difgk=0.0
	globalVar.tradat['difgk']=difgk 

	"""fp-concentration in graphite grains"""
	difkcp=[]
	idifkcp=5
	for i in range(idifkcp):
		difkcp.append(0.0)
	globalVar.tradat['difkcp']=difkcp

	"""fraction of defective particles"""
	pbra=0.0
	globalVar.cpbru['pbra']=pbra 

	"""number of intact particles"""
	pzahli=0
	globalVar.cpbru['pzahli']=pzahli

	"""number of defective particles"""
	pzahld=[]
	ipzahld=10
	for i in range(ipzahld):
		pzahld.append(0.0)
	globalVar.cpbru['pzahld']=pzahld

	"""fractional release from CP-kernel (avg. over all cp incl. defectives) """
	freik=0
	globalVar.cpbru['freik']=freik

	"""fractional release from particles (avg. over all cp incl. defectives) """
	frcp=0
	globalVar.cpbru['frcp']=frcp

	"""fractional release from fuel element"""
	frb=0
	globalVar.cpbru['frb']=frb

	"""fractional release from graphite grain"""
	frgk=0
	globalVar.cpbru['frgk']=frgk

	return
