import math,sys,shutil,os
from InputData import readdata

if __name__=="__main__":
	if len(sys.argv)==2:
		f=sys.argv[1]
	else:
		print("Argumen salah...")
		exit(1)
		
	data=readdata(f)
	print('ngnr='+str(data['ngnr']))
	print('title='+data['title'])
	print('te='+str(data['te']))
	print('dt='+str(data['dt']))
	print('dtout='+str(data['dtout']))
	print('uk='+str(data['uk']))
	print('rbe='+str(data['rbe']))
	print('nrbe='+str(data['nrbe']))
	print('pzahl0='+str(data['pzahl0']))
	print('zeit0d='+str(data['zeit0d']))
	print('xn0='+str(data['xn0']))
	print('zeitpr='+str(data['zeitpr']))
	print('gamma='+str(data['gamma']))
	print('rcp='+str(data['rcp']))
	print('rkorn='+str(data['rkorn']))
	print('nrc='+str(data['nrc']))
	print('nk0='+str(data['nk0']))
	print('unkongk='+str(data['unkongk']))
	print('unkontp='+str(data['unkontp']))
