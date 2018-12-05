import sys, math, re
def readdata(namafile):
	f=open(namafile,"r")
	data={}
	
	ngnr=int(f.readline())
	data['ngnr']=ngnr
	
	title=''
	for i in range(4):
		a=f.readline()
		title=title+(a.split('\n'))[0]
	a=title.split()
	title=''
	for i in range(len(a)):
		title=title+a[i]+' '
	data['title']=title
		
	te=[]
	temp1=f.readline()
	temp2=temp1.split()
	temp3=float(temp2[0])
	te.append(temp3)
	temp3=float(temp2[1])
	te.append(temp3)
	data['te']=te
	
	dt=[]
	temp1=f.readline()
	temp2=temp1.split()
	temp3=float(temp2[0])
	dt.append(temp3)
	temp3=float(temp2[1])
	dt.append(temp3)
	data['dt']=dt
	
	dtout=[]
	temp1=f.readline()
	temp2=temp1.split()
	temp3=float(temp2[0])
	dtout.append(temp3)
	temp3=temp2[1]
	dtout.append(temp3)
	data['dtout']=dtout
	
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
	data['uk']=uk
	
	rbe=[]
	temp1=f.readline()
	temp2=temp1.split()
	rbe.append(float(temp2[0]))
	rbe.append(float(temp2[1]))
	data['rbe']=rbe
	
	nrbe=[]
	temp1=f.readline()
	temp2=temp1.split()
	nrbe.append(int(temp2[0]))
	nrbe.append(int(temp2[1]))
	data['nrbe']=nrbe
	
	temp1=f.readline()
	temp2=temp1.split()
	data['pzahl0']=float(temp2[0])
	data['zeit0d']=float(temp2[1])
	data['xn0']=0
	data['zeitpr']=float(temp2[2])
	data['gamma']=float(temp2[3])
	
	rcp=[]
	temp1=f.readline()
	temp2=temp1.split()
	temp3=float(temp2[0])
	rcp.append(temp3)
	temp3=float(temp2[1])
	rcp.append(temp3)
	temp3=float(temp2[2])
	rcp.append(temp3)
	temp3=float(temp2[3])
	rcp.append(temp3)
	temp3=float(temp2[4])
	rcp.append(temp3)
	data['rcp']=rcp
	
	temp1=f.readline()
	temp2=temp1.split()
	data['rkorn']=float(temp2[0])
	
	nrc=[]
	temp1=f.readline()
	temp2=temp1.split()
	temp3=int(temp2[0])
	nrc.append(temp3)
	temp3=int(temp2[1])
	nrc.append(temp3)
	temp3=int(temp2[2])
	nrc.append(temp3)
	temp3=int(temp2[3])
	nrc.append(temp3)
	temp3=int(temp2[4])
	nrc.append(temp3)
	data['nrc']=nrc
	
	temp1=f.readline()
	temp2=temp1.split()
	data['nk0']=int(temp2[0])
	
	unkontp=[]
	temp1=f.readline()
	temp2=temp1.split()
	data['unkongk']=float(temp2[0])
	unkontp.append(float(temp2[1]))
	unkontp.append(float(temp2[2]))
	unkontp.append(float(temp2[3]))
	unkontp.append(float(temp2[4]))
	data['unkontp']=unkontp
	
	temp1=f.readline()
	temp2=temp1.split()
	temp3=len(temp2)
	data['d0g']=float(temp2[0])
	data['AKG']=float(temp2[1])
	data['F0G']=0.0
	data['F0PK']=0.0
	data['F0PSIC']=0.0
	data['F0PDPK']=0.0
		
		
	temp1=f.readline()
	temp2=temp1.split()
	data['d0k']=float(temp2[0])
	data['AKK']=float(temp2[1])
	data['F0KG']=0.0
	
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
	data['d0p']=d0p
	
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
	data['akp']=akp
	
	temp1=f.readline()
	temp2=temp1.split()
	data['recker']=float(temp2[0])
	data['recpyc']=float(temp2[1])
	data['recgra']=float(temp2[2])
	
	temp1=f.readline()
	temp2=temp1.split()
	data['aci']=float(temp2[0])
	data['bci']=float(temp2[1])
	data['cgrenc']=float(temp2[2])
	data['ifadc']=float(temp2[3])
	data['tgrenc']=0.0
	
	temp1=f.readline()
	temp2=temp1.split()
	data['anci']=float(temp2[0])
	data['bnci']=float(temp2[1])
	data['enci']=float(temp2[2])
	data['fnci']=float(temp2[3])
	
	temp1=f.readline()
	temp2=temp1.split()
	data['ainv']=float(temp2[0])
	if data['ngnr']==1:
		data['zerfk']=7.29E-10
	elif data['ngnr']==2:
		data['zerfk']=7.712E-10
	elif data['ngnr']==3:
		data['zerfk']=3.21E-8
	elif data['ngnr']==4:
		data['zerfk']=1.0E-6
	
	f.close()
	return data
