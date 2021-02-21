
import pandas as pd
import numpy as np
# convert to magnitudes
import numpy as np
from pyphot import LickLibrary
import pandas as pd
from astropy.io import fits
from matplotlib import pyplot as plt

#uses pyphot library
#Measures lick indexes with errors in dependence of radial velocity-- radial velocity must be measured first.
#you may need more than 2 radial mesurements, or more than 2 apertures in fits file it doesn't work with a single apperture.
#computes errors from cardiel et al 1998
#the case of logaritmic wavelenght is on, fitsfile must be in logaritmic scale and containing more than 2 apertures.
#Don't ask me how but it works.

def medirindices(fits_file):
	def lick_measure(fits_file,indice):
		def lickerrcoef(indice): #computes errors from cardiel et al 1998
			df=pd.read_csv('index_definition.csv')
			df.set_index("index", inplace=True)
			a1=(1/(df.loc[[indice],['dc']].values))
			a2=df.loc[[indice],['Lr']].values-df.loc[[indice],['Lc']].values
			a3=df.loc[[indice],['Lr']].values-df.loc[[indice],['Lb']].values
			A4=np.power((a2/a3),2)
			a5=1/df.loc[[indice],['db']].values
			a6=df.loc[[indice],['Lc']].values-df.loc[[indice],['Lb']].values
			a7=df.loc[[indice],['Lr']].values-df.loc[[indice],['Lb']].values
			A8=np.power(a6/a7,2)
			a9=1/df.loc[[indice],['dr']].values
			global c2
			global c1
			global c3
			c2=np.sqrt((a1)+(A4)*(a5)+(A8)*(a9)) #equation 44
			c1=(df.loc[[indice],['dc']].values)*c2 #equation 43
			c3=2.5*c2*np.log10(np.e)
		 	return c1,c2,c3

		R=np.loadtxt('R',unpack=True,usecols=[0]) #file of Radius points
		df=pd.read_csv('index_definition.csv') #needs this file with definitions
		df.set_index("index", inplace=True)
		lib = LickLibrary()
		df1=pd.DataFrame()
	#	files=raw_input('Name of the spectras to get magnitudes')
		hdu=fits.open(fits_file)       
		crval = hdu[0].header["CRVAL1"]
		naxis = hdu[0].header["NAXIS1"]
		crpix = hdu[0].header.get("CRPIX1", 0)
		cdelt = hdu[0].header["CDELT1"]
		ltv = hdu[0].header.get("LTV1", 0)
		f = hdu[0].data
		k = \
		crval + (np.arange(naxis) - crpix) * cdelt - ltv * cdelt
		vel=np.loadtxt('Velocities',unpack=True, usecols=[0]) #file of radial velocities
		c = 299792.458                  # speed of light in km/s
		z_central=vel[11]/c #central z
		zzz=[]
		df1['R']=R
		df1['vel']=vel
		df1['z']=vel/c
		if indice == 'Hb':
			indice2='H_beta'
		elif indice == 'HdA':
			indice2='Hdelta_A'
		elif indice == 'HdF':
			indice2='Hdelta_F'
		elif indice == 'HgA':
			indice2='Hgamma_A'
		elif indice == 'HgF':
			indice2='Hgamma_F'
		elif indice == 'NaD':
			indice2='Na_D'
		elif indice == 'Mgb':
			indice2='Mg_b'
		elif indice == 'Mg2':
			indice2='Mg_2'
		elif indice == 'Mg1':
			indice2='Mg_1'
		elif indice == 'CN1':
			indice2='CN_1'
		elif indice == 'CN2':
			indice2='CN_2'
		elif indice == 'C4668':
			indice2='Fe4668'
		else:
			indice2=indice

		zvar=[]
		men_flux=[]
		men_std=[]
		for i in range(0,len(f)):
			k =crval + (np.arange(naxis) - crpix) * cdelt - ltv * cdelt
			uu=10**k #case of logwavelenght
			#t=k/(1+z)
			uu=uu/(1+df1['z'][i])
			zvar.append(uu)
		w_reound=np.round(zvar)
		lickerrcoef(indice)
		for i in range(0,len(f)):
			try:
				men_flux.append(np.mean(f[i][np.argwhere((w_reound[i] == np.round(df.loc[[indice],['cb1']].values)))[0][1]:np.argwhere((w_reound[i] == np.round(df.loc[[indice],['cb2']].values)))[0][1]]))
				pass
			except Exception:
				men_flux.append(np.nan)
			try:
				men_std.append(np.std(f[i][np.argwhere((w_reound[i] == np.round(df.loc[[indice],['cb1']].values)))[0][1]:np.argwhere((w_reound[i] == np.round(df.loc[[indice],['cb2']].values)))[0][1]]))
				pass
			except Exception:
				men_std.append(np.nan)

		men_flux=np.array(men_flux)
		men_std=np.array(men_std)
		global snr
		global sigma_molecular
		snr=men_flux/men_std
		sigma_molecular=c3/snr #equation 45
		mm=lib[indice2]
		indexxX=[]
		snrrr=[]
		for i in range(0,len(f)):
		  index2=mm.get(zvar[i],f[i],axis=1)
		  indexxX.append(index2)
		global indexXxx
		indexXxx=np.array(indexxX)
		global sigma_atomic
		sigma_atomic=(c1-c2*indexxX)/snr
		#plt.errorbar(df1['R'],indexxX,yerr=sigma_atomic[0])
		#plt.show()


	indice=['Mg2','HgF','HgA','HdA','HdF','Ca4227','G4300','Fe4383','Hb','Mg1','Fe5015','Mgb','Fe5270','Fe5406','NaD','C4668','Fe5782','CN1','CN2','Fe5335']
	df_final=pd.DataFrame()
	R=np.loadtxt('R',unpack=True,usecols=[0])
	df_final['R']=R
	for i in range(0,len(indice)):
		lick_measure(fits_file,indice[i])
		df_final[indice[i]]=indexXxx

		if indice[i] == 'Mg1':
			df_final[indice[i]+'_err']=sigma_molecular[0]
		elif indice[i] == 'Mg2':
			df_final[indice[i]+'_err']=sigma_molecular[0]
		elif indice[i] == 'CN1':
			df_final[indice[i]+'_err']=sigma_molecular[0]
		elif indice[i] == 'CN2':
			df_final[indice[i]+'_err']=sigma_molecular[0]
		else:
			df_final[indice[i]+'_err']=sigma_atomic[0]
		df_final[indice[i]+'_snr']=snr
	df_final['<Fe>']=(( df_final['Fe5335']+ df_final['Fe5270'])/2)
	df_final['<Fe>_err']=np.sqrt((np.power( df_final['Fe5335_err'].astype(float),2)+np.power( df_final['Fe5270_err'].astype(float),2))/2)
	df_final['<Fe>_snr']=( df_final['Fe5335_snr'].astype(float)+ df_final['Fe5270_snr'].astype(float)/2)
	df_final['MgFe']=np.sqrt( df_final['<Fe>']* df_final['Mgb'])
	df_final['MgFe_err']=(( df_final['Mgb_err']/(2* df_final['Mgb'])))+(( df_final['<Fe>_err']/(2* df_final['<Fe>']))).astype(float)
	df_final['MgFe_snr']=( df_final['Fe5335_snr'].astype(float)+ df_final['Fe5270_snr'].astype(float)+ df_final['Mgb_snr'].astype(float))/3
	# df_final['[MgFe]']=mgfe
	df_final['Log(Mgb/<Fe>)']=np.log10( df_final['Mgb']/ df_final['<Fe>'])
	df_final['Log(Mgb/<Fe>)_err']=abs( df_final['Log(Mgb/<Fe>)'])*np.sqrt(( df_final['Mgb_err']/abs( df_final['Mgb']))**2+( df_final['<Fe>_err']/abs(df_final['<Fe>']))**2)
	name=fits_file.replace('.fits','')
	df_final.to_csv(str(name)+'calibred_index.csv',sep=',',index=True)





def medirlick(fits_file):
	def lick_measure(fits_file,indice):
		def lickerrcoef(indice):
			df=pd.read_csv('index_definition.csv')
			df.set_index("index", inplace=True)
			a1=(1/(df.loc[[indice],['dc']].values))
			a2=df.loc[[indice],['Lr']].values-df.loc[[indice],['Lc']].values
			a3=df.loc[[indice],['Lr']].values-df.loc[[indice],['Lb']].values
			A4=np.power((a2/a3),2)
			a5=1/df.loc[[indice],['db']].values
			a6=df.loc[[indice],['Lc']].values-df.loc[[indice],['Lb']].values
			a7=df.loc[[indice],['Lr']].values-df.loc[[indice],['Lb']].values
			A8=np.power(a6/a7,2)
			a9=1/df.loc[[indice],['dr']].values
			global c2
			global c1
			global c3
			c2=np.sqrt((a1)+(A4)*(a5)+(A8)*(a9)) #equation 44
			c1=(df.loc[[indice],['dc']].values)*c2 #equation 43
			c3=2.5*c2*np.log10(np.e)
		 	return c1,c2,c3


		df=pd.read_csv('index_definition.csv')
		df.set_index("index", inplace=True)
		lib = LickLibrary()
		df1=pd.DataFrame()
	#	files=raw_input('Name of the spectras to get magnitudes')
		hdu=fits.open(fits_file)       
		crval = hdu[0].header["CRVAL1"]
		naxis = hdu[0].header["NAXIS1"]
		crpix = hdu[0].header.get("CRPIX1", 0)
		cdelt = hdu[0].header["CDELT1"]
		ltv = hdu[0].header.get("LTV1", 0)
		f = hdu[0].data
		k = \
		crval + (np.arange(naxis) - crpix) * cdelt - ltv * cdelt
		vel,disp,h3,h4,h5,h6=np.loadtxt('Resultados/'+fits_file+'_Velocities_out',unpack=True, usecols=[0,1,2,3,4,5])
		c = 299792.458                  # speed of light in km/s

		if indice == 'Hb':
			indice2='H_beta'
		elif indice == 'HdA':
			indice2='Hdelta_A'
		elif indice == 'HdF':
			indice2='Hdelta_F'
		elif indice == 'HgA':
			indice2='Hgamma_A'
		elif indice == 'HgF':
			indice2='Hgamma_F'
		elif indice == 'NaD':
			indice2='Na_D'
		elif indice == 'Mgb':
			indice2='Mg_b'
		elif indice == 'Mg2':
			indice2='Mg_2'
		elif indice == 'Mg1':
			indice2='Mg_1'
		elif indice == 'CN1':
			indice2='CN_1'
		elif indice == 'CN2':
			indice2='CN_2'
		elif indice == 'C4668':
			indice2='Fe4668'
		else:
			indice2=indice

		men_flux=[]
		men_std=[]

		k =crval + (np.arange(naxis) - crpix) * cdelt - ltv * cdelt
		uu=10**k #case of logwavelenght
			#t=k/(1+z)
		uu=uu/(1+(vel/c))
		zvar=uu
		w_reound=np.round(zvar)
		lickerrcoef(indice)
#		for i in range(0,len(w_reound)):
		try:
			men_flux.append(np.mean(f[np.argwhere((w_reound == np.round(df.loc[[indice],['cb1']].values)))[0][1]:np.argwhere((w_reound == np.round(df.loc[[indice],['cb2']].values)))[0][1]]))
			pass
		except Exception:
			men_flux.append(np.nan)
		try:
			men_std.append(np.std(f[np.argwhere((w_reound == np.round(df.loc[[indice],['cb1']].values)))[0][1]:np.argwhere((w_reound == np.round(df.loc[[indice],['cb2']].values)))[0][1]]))
			pass
		except Exception:
			men_std.append(np.nan)

		men_flux=np.array(men_flux)
		men_std=np.array(men_std)
		global snr
		global sigma_molecular
		snr=men_flux/men_std
		sigma_molecular=c3/snr #equation 45
		mm=lib[indice2]
		#indexxX=[]
		snrrr=[]
		index2=mm.get(zvar,f,axis=1)
		indexxX=index2#.append(index2)
		global indexXxx
		indexXxx=np.array(indexxX)
		global sigma_atomic
		sigma_atomic=(c1-c2*indexxX)/snr
		#plt.errorbar(df1['R'],indexxX,yerr=sigma_atomic[0])
		#plt.show()


	indice=['Mg2','HgF','HgA','HdA','HdF','Ca4227','G4300','Fe4383','Hb','Mg1','Fe5015','Mgb','Fe5270','Fe5406','NaD','C4668','Fe5782','CN1','CN2','Fe5335']
	df_final=pd.DataFrame()
	df_final['id']=1
	df_final['vel']=[vel]
	for i in range(0,len(indice)):
		lick_measure(fits_file,indice[i])
		df_final[indice[i]]=[indexXxx]

		if indice[i] == 'Mg1':
			df_final[indice[i]+'_err']=sigma_molecular[0]
		elif indice[i] == 'Mg2':
			df_final[indice[i]+'_err']=sigma_molecular[0]
		elif indice[i] == 'CN1':
			df_final[indice[i]+'_err']=sigma_molecular[0]
		elif indice[i] == 'CN2':
			df_final[indice[i]+'_err']=sigma_molecular[0]
		else:
			df_final[indice[i]+'_err']=sigma_atomic[0]
		df_final[indice[i]+'_snr']=snr
	df_final['<Fe>']=(( df_final['Fe5335']+ df_final['Fe5270'])/2)
	df_final['<Fe>_err']=np.sqrt((np.power( df_final['Fe5335_err'].astype(float),2)+np.power( df_final['Fe5270_err'].astype(float),2))/2)
	df_final['<Fe>_snr']=( df_final['Fe5335_snr'].astype(float)+ df_final['Fe5270_snr'].astype(float)/2)
	df_final['MgFe']=np.sqrt( ((df_final['<Fe>']*df_final['Mgb']).values).astype(np.float64))
	df_final['MgFe_err']=(( df_final['Mgb_err']/(2* df_final['Mgb'])))+(( df_final['<Fe>_err']/(2* df_final['<Fe>']))).astype(float)
	df_final['MgFe_snr']=( df_final['Fe5335_snr'].astype(float)+ df_final['Fe5270_snr'].astype(float)+ df_final['Mgb_snr'].astype(float))/3
	# df_final['[MgFe]']=mgfe
	df_final['Log(Mgb/<Fe>)']=np.log10( (df_final['Mgb'].values).astype(np.float64)/ (df_final['<Fe>'].values).astype(np.float64))
	df_final['Log(Mgb/<Fe>)_err']=abs( df_final['Log(Mgb/<Fe>)'])*np.sqrt(( (df_final['Mgb_err'].values).astype(np.float64)/abs( (df_final['Mgb'].values).astype(np.float64)))**2+( (df_final['<Fe>_err'].values).astype(np.float64)/abs((df_final['<Fe>'].values).astype(np.float64)))**2)
	name=fits_file.replace('.fits','')
	df_final.to_csv(str(name)+'_indexes.csv',sep=',',index=True)
