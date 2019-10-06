import sys, math, re
def readdata(namafile):
    f=open(namafile,"r")
    global data
    data={}
    diff={}     # data terkait transport
    kont={}     # data terkait Heavy Metal Contamination
    recr={}     # data terkait recoil
    nukdat={}   # data terkait Data Nuclear
    adsoc={}    # data terkait Sorption Isotherms
    geod={}     # data terkait geometri
    geodi={}    # data terkait geometri
    sic={}      # data terkait Lapisan Silicon Carbide
    zeits={}    # data terkait waktu
    real1={}    # data terkait waktu dan fast flux
    itemp={}    # data terkait info jumlah data pd kondisi operasi normal dan kecelakaan
    ttemp={}    # detail data history temp pada kondisi operasi normal dan kecelakaan
    invent={}   # data terkait fp-inventory
    inven1={}   # data terkait fp-inventory
    konz={}     # data terkait konsentrasi FP
    tradat={}
    cpbru={}
    #
    #
    # Start Reading Input File
    #
    ngnr=int(f.readline())
    data['ngnr']=ngnr
	
    title=''
    for i in range(4):
        title=title+f.readline()
    data['title']=title
		
    te=[]
    temp1=f.readline()
    temp2=temp1.split()
    item=len(temp2)
    for i in range(item):
        te.append(float(temp2[i]))
    zeits['te']=te
	
    dt=[]
    temp1=f.readline()
    temp2=temp1.split()
    for i in range(item):
        dt.append(float(temp2[i]))
    zeits['dt']=dt
	
    dtout=[]
    temp1=f.readline()
    temp2=temp1.split()
    item=len(temp2)
    for i in range(item):
        dtout.append(float(temp2[i]))
    zeits['dtout']=dtout
    #
    data['zeits']=zeits
    #
	
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
    item=len(temp2)
    for i in range(item):
        rbe.append(float(temp2[i]))
    geod['rbe']=rbe                     # outer radius of i-th zone in the in the fuel element (or pebble fuel).
    
	
    nrb=[]
    temp1=f.readline()
    temp2=temp1.split()
    item=len(temp2)
    npb=0
    for i in range(item):
        nrb.append(int(temp2[i]))
        npb=npb+int(temp2[i])
    geodi['nrb']=nrb                    # number of mesh in each zone
    geodi['nrbe']=item                  # number of total zones in the fuel element
    geodi['npb']=npb                    # number of total meshes in fuel element (adding meshes from all zones)
    
    temp1=f.readline()
    temp2=temp1.split()
    geod['pzahl0']=float(temp2[0])      #number of particles in fuel sphere
    real1['zeit0d']=float(temp2[1])      #normal operation/irradiation time
    real1['xn0']=0                       #total inventory of the fission product species considered
    real1['zeitpr']=float(temp2[2])      #time after begin of the accident/heating phase.
    real1['gamma']=float(temp2[3])       # fast neutron fluence at begin of the accident/heating phase
    #
    data['real1']=real1
    #
    rcp=[]                          
    temp1=f.readline()              # first rcp, rcp[0] is the outer diam. of fuel kernel.
    temp2=temp1.split()
    item=len(temp2)
    for i in range(item):
        rcp.append(float(temp2[i]))
    geod['rcp']=rcp                 # outer radius of i-th zone in the in the coated particle.
	
    temp1=f.readline()
    temp2=temp1.split()
    geod['rkorn']=float(temp2[0])       # radius of grahite grain (graphite pore is the grain boundary)
	
    nrc=[]
    temp1=f.readline()
    temp2=temp1.split()
    item=len(temp2)
    npp=0
    for i in range(item):
        nrc.append(int(temp2[i]))
        npp=npp+int(temp2[i])
    geodi['nrc']=nrc                # number of mesh in each zone
    geodi['nrp']=item               # number of zones in fuel particle
    geodi['npp']=npp                # number of total mesh in fuel particle (adding meshes in all zones)
     
    temp1=f.readline()
    temp2=temp1.split()
    geodi['nk0']=int(temp2[0])
	 #
    #Geometry data calculated from raw input
    #
    # zone radius in fuel element : array 'rb'
    rb=[]
    rb.append(0.0)                  # it's the center of the fuel element, as first rb array.
    for i in range(1,npb+1):        # warning: current algorithm is only suitable for nrbe=2
        if i<=nrb[0]:
            mesh=rbe[0]/nrb[0]
            #print(i, rb[i-1], mesh)
            rb.append(rb[i-1]+mesh)
        else:
            mesh=(rbe[1]-rbe[0])/nrb[1]
            #print(i, rb[i-1], mesh)
            rb.append(rb[i-1]+mesh)
    #print('-----')
    #print('radius of zone in fuel element: ', len(rb),rb)
    geod['rb']=rb
    #
    # zone radius in fuel particle : array 'rp'
    mesh=[]
    accitem=0
    acc=[]
    nrp=geodi['nrp']
    for i in range(nrp):
        if i==0:
            mesh.append(rcp[0]/nrc[0])
            accitem=accitem+nrc[0]
            acc.append(accitem)
        else:
            mesh.append((rcp[i]-rcp[i-1])/nrc[i])
            accitem=accitem+nrc[0]
            acc.append(accitem)
    #print('mesh in fuel particle: ', len(mesh),mesh)
    #print('accumulative nrc : ' , len(acc), acc)
    #  
    rp=[]
    rp.append(0.0)                  # it's the center of the fuel particle, as first rp array.
    temp=0.0
    for i in range(nrp):        # warning: current algorithm is only suitable for nrbe=2
        for j in range(nrc[i]):
            rp.append(temp+mesh[i])
            temp=temp+mesh[i]
    print('-----')
    #print('radius of zone in fuel particle: ')
    #for i in range(npp+1):
    #    print(i, rp[i])
    #print('-----')
    geod['rp']=rp
    #
    #
    # zone radius in graphite grain
    rk0=[]
    rk0.append(0.0)
    nk0=geodi['nk0']
    rkorn=geod['rkorn']
    for i in range(1,nk0+1):
        rk0.append(rk0[i-1]+(rkorn/nk0))
    #print('radius in graphite grain: ',len(rk0))
    #print('nk0, rkorn, rkorn/nk0 : ', nk0, rkorn, rkorn/nk0)
    #print(rk0)
    geod['rk0']=rk0
    #
    data['geod']=geod 
    data['geodi']=geodi 
    #
    #   Heavy Metal Contamination Data
    #
    ukontp=[]
    temp1=f.readline()
    temp2=temp1.split()
    kont['ukongk']=float(temp2[0])
    kont['ukongp']=float(temp2[1])
    unkernel=1.-(float(temp2[0])+float(temp2[1])+float(temp2[2])+float(temp2[3])+float(temp2[4])+float(temp2[5]))
    ukontp.append(unkernel)
    ukontp.append(float(temp2[2]))
    ukontp.append(float(temp2[3]))
    ukontp.append(float(temp2[4]))
    ukontp.append(float(temp2[5]))
    kont['ukontp']=ukontp
    #
    data['kont']=kont
    
    #
    #   Transport Data
    #
    temp1=f.readline()
    temp2=temp1.split()
    temp3=len(temp2)
    diff['d0g']=float(temp2[0])
    diff['akg']=float(temp2[1])
    diff['f0g']=float(temp2[2])
    diff['f0pk']=float(temp2[3])
    diff['f0psic']=float(temp2[4])
    diff['f0pdpk']=float(temp2[5])
		
    temp1=f.readline()
    temp2=temp1.split()
    diff['d0k']=float(temp2[0])
    diff['akk']=float(temp2[1])
    diff['f0kg']=float(temp2[2])
    
    	
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
    diff['d0p']=d0p
	
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
    diff['akp']=akp
    
    data['diff']=diff

    temp1=f.readline()
    temp2=temp1.split()
    recr['recker']=float(temp2[0])
    recr['recpyc']=float(temp2[1])
    recr['recgra']=float(temp2[2])
    #
    data['recr']=recr
    #	
    
    temp1=f.readline()
    temp2=temp1.split()
    adsoc['aci']=float(temp2[0])
    adsoc['bci']=float(temp2[1])
    adsoc['cgrenc']=float(temp2[2])
    adsoc['ifadc']=float(temp2[3])
    adsoc['tgrenc']=0.0
	
    temp1=f.readline()
    temp2=temp1.split()
    adsoc['anci']=float(temp2[0])
    adsoc['bnci']=float(temp2[1])
    adsoc['enci']=float(temp2[2])
    adsoc['fnci']=float(temp2[3])
    #
    data['adsoc']=adsoc
    #
    temp1=f.readline()
    temp2=temp1.split()
    nukdat['ainv']=float(temp2[0])
    if data['ngnr']==1:
        nukdat['zerfk']=7.29E-10
    elif data['ngnr']==2:
        nukdat['zerfk']=7.712E-10
    elif data['ngnr']==3:
        nukdat['zerfk']=3.21E-8
    elif data['ngnr']==4:
        nukdat['zerfk']=1.0E-6
    #
    data['nukdat']=nukdat
    #
    temp1=f.readline()      #to read the 'Irradiation Temp. History line
    #
    ####
    # Read Normal Operation / Irradiation Data ( time , temp)
    ####
    temp1=f.readline()
    temp2=temp1.split()
    ni0=int(temp2[0])               # ni0 : number of pair of data in normal operation
    itemp['ni0']=ni0
    #print('ni0 : ', ni0)
    azeit0=[]
    atemp0=[]
    for i in range(ni0):
        temp1=f.readline()
        temp2=temp1.split()
        azeit0.append(float(temp2[1]))
        atemp0.append(float(temp2[2]))
    ttemp['azeit0']=azeit0          # units in [days]
    ttemp['atemp0']=atemp0          # units in degC
    #print('azeit0 : ', len(azeit0),azeit0)
    #print('----')
    #print('atemp0 : ', len(atemp0),atemp0)
    #print('----')
    # data terkait Lapisan Silicon Carbide
    temp1=f.readline()      #to read the 'Accident Temp. History line
    #
    ####
    # Read Normal Operation / Irradiation Data ( time , temp)
    ####
    temp1=f.readline()
    temp2=temp1.split()
    ni=int(temp2[0])
    itemp['ni']=ni
    #print('ni0 : ', ni0)
    azeit=[]
    atemp=[]
    for i in range(ni):
        temp1=f.readline()
        temp2=temp1.split()
        azeit.append(float(temp2[1]))
        atemp.append(float(temp2[2]))
    ttemp['azeit']=azeit                #units in [hours]
    ttemp['atemp']=atemp                #units in degC
    #print('azeit : ', len(azeit),azeit)
    #print('----')
    #print('atemp : ', len(atemp),atemp)
    #print('----')
    #
    data['itemp']=itemp
    data['ttemp']=ttemp
    
    ####
    # data terkait Lapisan Silicon Carbide : SIC
    ####
    alsic=0                 #corroded length of SiC layer in cm
    ds0=rcp[3]              #thickness of SiC layer , taken from rcp input array.
    dkor=0.0                # ?
    isic=nrc[3]             # number of SiC-layer
    dsic=ds0/isic           # thickness of a spherical shell in the SiC-layer.
    dif=0.0
    msic=0                  #number of corroded sphere shells in SiC
    msiges=0                # 
    sic['alsic']=alsic
    sic['ds0']=ds0
    sic['dkor']=dkor
    sic['dsic']=dsic
    sic['dif']=dif
    sic['isic']=isic
    sic['msic']=msic
    sic['msiges']=msiges
    #
    data['sic']=sic
    #
    #
    ####
    # data terkait profile konsentrasi: INVENT
    ####
    #
    aicp=[]                     # fp-inventory in single layers 
    for i in range(5):
        aicp.append(0.0)
    invent['aicp']=aicp
    #
    aicpi=0.0 
    invent['aicpi']=aicpi       # fp-inventory in intact coated particle (CP)
    #
    aicpd=[]
    for i in range(5):
        aicpd.append(0.0)
    invent['aicpd']=aicpd       # sp-inventory in one defective CP (kernel)
    #
    aicpk=0.0
    invent['aicpk']=aicpk       # fp-inventory in kernels of all coated particle (CP)
    #
    aicpdg=0.0
    invent['aicpdg']=aicpdg     # fp-inventory in kernels of all defective coated particle (CP)
    #
    aicpg=0.0
    invent['aicpg']=aicpg       # fp-inventory in all coated particles (CPs)
    #
    aigk=0.0
    invent['aigk']=aigk         # fp-inventory in dern graphitekoernern
    #
    data['invent']=invent
    #
    ####
    # data terkait profile inventory: INVEN1
    ####
    #
    aigp=0.0
    inven1['aigp']=aigp         # fp-inventory in graphite-pores
    #
    aifr=0.0
    inven1['aifr']=aifr         # fp-amount outside the fuel element 
    #
    data['inven1']=inven1
    #
    ####
    # data terkait konsentrasi FP: KONZ
    ####
    cb=[]
    icb=200
    for i in range(icb):
        cb.append(0.0)
    konz['cb']=cb       # fp-concentration in graphite pores
    #
    ck=[]
    ick=200
    for i in range(ick):
        ck.append(0.0)
    konz['ck']=ck       # fp-concentration in graphite grains
    #
    cp=[]
    icp=200
    for i in range(icp):
        cp.append(0.0)
    konz['cp']=cp       # fp-concentration in graphite particles
    #
    cpb=[]
    icpb=200
    for i in range(icpb):
        cpb.append(0.0)
    konz['cpb']=cpb       # fp-concentration in graphite grains
    #
    data['konz']=konz
    #
    ####
    # data terkait ....: TRADAT
    ####
    temper=0.0
    tradat['temper']=temper         # temperature 
    #
    difgp=[]
    idifgp=5
    for i in range(idifgp):
        difgp.append(0.0)
    tradat['difgp']=difgp       # fp-concentration in graphite grains
    #
    difgk=0.0
    tradat['difgk']=difgk         # fp-amount outside the fuel element 
    #
    difkcp=[]
    idifkcp=5
    for i in range(idifkcp):
        difkcp.append(0.0)
    tradat['difkcp']=difkcp       # fp-concentration in graphite grains
    #
    data['tradat']=tradat
    #
    ####
    # data terkait ....: CPBRU
    ####
    #
    pbra=0.0
    cpbru['pbra']=pbra         # fraction of defective particles 
    #
    pzahli=0
    cpbru['pzahli']=pzahli         # number of intact particles 
    #
    pzahld=[]
    ipzahld=10
    for i in range(ipzahld):
        pzahld.append(0.0)
    cpbru['pzahld']=pzahld       # number of defective particles
    #
    freik=0
    cpbru['freik']=freik         # fractional release from CP-kernel (avg. over all cp incl. defectives) 
    #
    frcp=0
    cpbru['frcp']=frcp         # fractional release from particles (avg. over all cp incl. defectives) 
    #
    frb=0
    cpbru['frb']=frb         # fractional release from fuel element
    #
    frgk=0
    cpbru['frgk']=frgk         # fractional release from graphite grain
    #
    data['cpbru']=cpbru
    ##
    f.close()
    return data