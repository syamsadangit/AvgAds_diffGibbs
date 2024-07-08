import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random as rd
import re 
from sys import exit
import pickle
import warnings

warnings.filterwarnings("default", ".*", FutureWarning, r".*", 0)


pd.options.mode.chained_assignment = None

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#DEFINE FUNCTIONS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def iseven(num):
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#CHECKS IF num IS EVER OR NOT AND RETURNS condition
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if (num % 2 == 0):
		condition = True
	else:
		condition = False
	return condition

def FreeEntVib(T,Evib):
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#INPUT TEMPERATURE T AND VIBRATIONAL ENERGY EVIB
	#OUTPUT FREE ENERGY OF VIBRATION AND VIB TERM, ENT TERM
	#FOR REF https://doi.org/10.1103/PhysRevB.65.035406
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	kb = 8.617333262e-05 #eV/K 
	B = 1/(kb*T)
	x = B*Evib
	exp = np.exp(2*x)
	termvib     = Evib + (2*Evib/(exp - 1))
	termentropy = (1/B)*( (2*x/(exp - 1)) - np.log(1-(1/exp)) )   
	F = termvib - termentropy
	return F,termvib,termentropy

def write_vib_entropy(Vibenefile):
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#WRITING VIB FREE ENERGY AND ENT DUE TO VIB TO FILE
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	head_file = "Entropy_correction"
	dictinp = {}
	dfvib,dictvib = VibEneCollect(Vibenefile,dictinp)
	col = dfvib.columns
	Evib = dfvib[col[0]].copy() #SettingWithCopyWarning
	Free_energy_entropy,termvib,termentropy = FreeEntVib(T,Evib)
	dfvib.insert(1,'Entrene',termentropy,True)
	dfvib.insert(2,'tot_corr',Free_energy_entropy,True)
	fname  = Vibenefile.split(".")
	svname = head_file+"".join(fname[:-1])+".dat"
	dfvib.to_csv(svname,sep='\t',index=False,\
	float_format='%5.8f')
	return
	 
			
		

def EneCorrectVib(Vibenefile,dfEads,Eslab):
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#THIS FUNC TAKES ENERGY OF THE OPTIMISED ADSORBED SLAB
	#ENERGY OF REFERENCE SLAB, VIB ENERGY FILE 
	#AND ADD THE CORRESPONDING SITE-VIBRATIONAL ENERGY
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	print("The vib ene data is taken from"+Vibenefile)
	PHlst = ["PH2","PH3"] # UPDATE THIS TO MORE CONFIG IF NEEDED
	dictinp = {}
	dfvib,dictvib = VibEneCollect(Vibenefile,dictinp) #CALL FUNCTION
	remkeys = dictvib.keys()
	for ky in remkeys:
		Evib = dictvib[ky]
		dictvib[ky],V,E =  FreeEntVib(T,Evib) #CALL FUNCTION

	if Eslab == EslabA7P:
		Eslab = Eslab+dictvib["Ni3P2+4P"]
	elif Eslab == EslabAA7:
		Eslab = Eslab
	else:
		print("Check the slab energy and vibrational energies") 
	dfeads = dfEads.copy()
	for i in range(len(dfeads.iloc[:,2])):
		rem = dfeads.iloc[i,2] #col 2 is Remarks
		print("__________________________________________________")
		print("		{} row Remarks = ".format(i),rem)
		config = re.split("#|-",rem)
		for val in config:
			if len(val)==0 or val.isdigit():
				config.remove(val)  
		print("split remarks with H = ",config)
		for Helem in config:
			nH_con = re.split("@",Helem)
			print("split remarks = ",nH_con)
			gen_rem   = nH_con[-1]
			if gen_rem in PHlst:
				nPH = re.split("H",gen_rem)	
				div = int(nPH[-1])
				mltpl = int(nH_con[0])/div
			else:
				mltpl = int(nH_con[0])
			elem = nH_con[1] 
			if elem in dictvib:
				viben = mltpl*dictvib[elem]
				print("There are {} no \
				of H at the site ".format(mltpl)\
				,elem)
				print("Vib ene to add = {}".format(viben))	
				to_rep = dfeads.iloc[i,1].copy()
				print("Bef corr Eads = {}".format(to_rep))
				rep    = to_rep + viben
				print("Aft corr Eads = {}".format(rep))
				dfeads.at[i,"Eads"]=rep
				print("Corrected Value = {}".format(dfeads.iloc[i,1]))
			else:
				print("The configuration "+elem+" has no vibrational energy in the file "+Vibenefile)
		print("--------------------------------------------------")
	return dfeads,Eslab


def VibEneCollect(input_file,dict):
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#READS IN INPUT FILE CONTAINIG VIB ENERGY AS DF AND DICT 
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	col1 = "VibEne"
	col2 = "Remarks"
	df = pd.read_csv(input_file,skiprows=1,sep=" ",\
		names=(col1,col2))
	for i in range(len(df[col2])):
		key  = df[col2][i].split("#")[1]
		dict["{}".format(key)] = df[col1][i]
	print(" ")
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
	print("Vibrational Energy Colleced from {}".format(input_file))
	print(dict)
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
	print(" ")
	return df,dict


def AvgEads(df,Eslab,Vibenefile):
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#CALCULATES ADS ENERGY PER H AND GIBBS FREE ENERGY PER H
	#ADS_ENE_PER_H  = (E_slab+H - E_slab - (1/2 * n * EH2) )
	#DF1 IS ADS ENERGY PER H (NO CORRECTION)
	#DF2 IS GIBBS FREE ENERGY PER H (VIB-TS CORRECTED)
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	df1 = df.copy()
	y = (df1.Eads - Eslab - (df1.nH*EH2*0.5))/df1.nH
	dfeads = df.copy()
	df2,Eslab = EneCorrectVib(Vibenefile,dfeads,Eslab) #CALL FUNCTION
	y2 = (df2.Eads - Eslab - (df2.nH*EH2tot*0.5))/df2.nH
	df1.insert(2,'Eadsavg',y,True)
	df2.insert(2,'Eadsavg',y2,True)
	return df1,df2 

def diffGibbs(df1,df2,neibr):
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	#CALCULATES THE DIFFERENTIAL GIBBS FREE ENERGY
	#USE DF1,DF2 = CLEANDF(DF) TO OBTAIN DF1 AND DF2
	#DF1 & DF2 ARE ADS ENERGY VALUES AFTER CORRECTIONS
	#DF1 - WITHOUT REF ENE VAL. DF2 - WITH REF ENE VAL 
	#DELGN=EN - EN-1 - 0.5*EH2TOT + E*U (U=0 HERE)
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dropind = len(df2.index) - 1
	dfn     = df1.copy()
	dfn_1   = df2.copy()
	arn     = dfn.values.tolist()   
	arn_1   = dfn_1.values.tolist() 
	terms   = (-0.5*EH2tot*neibr) + eU
	dgibbs  = []
	for i in arn:
		for j in arn_1:
			nHin   = i[0]
			nHjn_1 = j[0]
			if (neibr==2):
				cond1=iseven(nHin)
				cond2=iseven(nHjn_1)
			elif (neibr==1):
				cond1 = True
				cond2 = True
			else:
				exit("Code not updated for neibr values other than 1 or 2")	
			Ein    = i[1]
			Ejn_1  = j[1]
			remin  = i[2]
			remjn_1= j[2]
			print(cond1, cond2)
			print(nHin,nHjn_1,nHin-nHjn_1,neibr)
			if(cond1 and cond2 and (nHin-nHjn_1)==neibr):
				dE   = Ein - Ejn_1 + terms
				rem  = remin+","+remjn_1 
				updt = [nHin,Ein,Ejn_1,dE,rem]
				dgibbs.append(updt)	
	DGibbsdf  = pd.DataFrame(dgibbs,columns=["nH","En",\
	"En-1","diffGibbs","Remarks(n,n-1)"])
	return DGibbsdf 
	 

def cleandf(df,Eslab,Vibenefile,rem_avail=False):
	# Finds the minimmum energy adsorption configuration
	# for a given number of Hydrogen. Then it discards 
	# all the configuration except the minimum energy 
	#configuration
	E_RT = 0 #0.03  
	df3 = df.copy()
	df3,Eslab = EneCorrectVib(Vibenefile,df3,Eslab)
	if(rem_avail):
		#REMG IS A GLOBAL VARIABLE DEFINED IN FUNCTION plotAvgAdsene
		drpind = df3.loc[ ~df3["Remarks"].isin(remG) ].index	
		df3.drop(drpind,inplace=True)
	dfsort = df3.sort_values(by=["nH","Eads"])
	#dupbool = dfsort.duplicated(subset=["nH"])
	#dfdrop = dfsort.drop_duplicates(subset=["nH"],\
	#keep='first',ignore_index=True)
	#print(dfdrop)
	dfdrop = dfsort.copy()
	nHEads0  = dfdrop[["nH","Eads","Remarks"]].values.tolist()
	nHEads1  = dfsort[["nH","Eads","Remarks"]].values.tolist()
	nHEads2  = []
	for i in nHEads0:
		for j in nHEads1:
			nHi0 = i[0]
			nHj1 = j[0]
			Ei0  = i[1]
			Ej1  = j[1]
			remi0= i[2]
			remj1= j[2]
			if(nHi0==nHj1):
				if(abs(Ej1-Ei0)<=E_RT):
					updt = [nHj1,Ej1,remj1]
					nHEads2.append(updt)
	dfdr1 = pd.DataFrame(nHEads2,columns=["nH","Eads","Remarks"])
	print("cleandf")
	print(dfsort)
	print(dfdr1)
	col = df3.columns
	dfref = pd.DataFrame([[0 ,Eslab,"#Refene"]],\
	columns = col,index = ["ref"] )
	dfclean = dfref.append(dfdr1,ignore_index=True)
	return dfdr1, dfclean

def write_AvgEads(input_file,Eslab,system,Vibenefile):
	df = pd.read_csv(input_file,skiprows=1,sep=" ",\
	names=('nH','Eads','Remarks')) 
	df1,df2 = AvgEads(df,Eslab,Vibenefile)   #CALL FUNCTION
	plotAvgAdsene(df1,df2,input_file,system) #CALL FUNCTION
	fname  = input_file.split(".")
	svname1 = "".join(fname[:-1])+saveAvg
	df1.to_csv(svname1,sep='\t',index=False,\
	float_format='%5.8f')
	svname2 = "".join(fname[:-1])+saveGibbs
	df2.to_csv(svname2,sep='\t',index=False,\
	float_format='%5.8f')
	print("\n")
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
	print("Input file "+input_file+" read.")
	print("Average Adsorption Energy \
	written to ---> "+svname1)
	print("Average Adsorption Energy \
	written to ---> "+svname2)
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
	print("\n")
	return 

def write_diffGibbs(input_file,Eslab,system,Vibenefile):
	df = pd.read_csv(input_file,skiprows=1,sep=" ",\
	names=('nH','Eads','Remarks'))
	dfcp = df.copy()
	dfpermucp = df.copy()
	df2,df3 = cleandf(dfcp,Eslab,Vibenefile,True)
	df2permu,df3permu=cleandf(dfcp,Eslab,Vibenefile,False)
	df4 = diffGibbs(df2,df3,1)
	df5 = diffGibbs(df2permu,df3permu,1)
	plotdGibbs(df4,input_file,system)
	plotdGibbs(df5,"permu"+input_file,system)
	plotdAdsEne(df4,input_file,system)
	fname  = input_file.split(".")
	svname = "".join(fname[:-1])+savedGib
	df4.to_csv(svname,sep='\t',index=False,\
	float_format='%5.8f')
	df5.to_csv("permu"+svname,sep='\t',index=False,\
	float_format='%5.8f')
	print("\n")
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
	print("Input file "+input_file+" read.")
	print("Differential Gibbs Adsorption Energy \
	written to ---> "+svname)
	print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
	print("\n")
	return

def plotAvgAdsene(dfEads,df,input_file,system):
	global remG
	fname  = input_file.split(".")
	minsvname = "".join(fname[:-1])+".dat"
	svname = "".join(fname[:-1])+pltEavg 
	pos = 50
	choice = np.arange(1,pos,1) 
	x = df['nH']
	y = dfEads['Eadsavg']#*df.nH
	#vibTScorre = dZPE - TScorr
	#y1 = y - TScorr + eU
	y1= df['Eadsavg'] + eU
	fig, (ax1,ax2) = plt.subplots(1,2,figsize=(20, 15))
	df1 = df.sort_values(by=["nH","Eadsavg"])
	df2 = df1.drop_duplicates(subset=["nH"],\
        keep='first',ignore_index=True)
	df11 = dfEads.sort_values(by=["nH","Eadsavg"])
	df21 = df11.drop_duplicates(subset=["nH"],\
        keep='first',ignore_index=True)
	#print(df2)
	xmin,ymin,rem = df21["nH"], df21["Eadsavg"] ,df21["Remarks"]
	df21.to_csv("AdsFenemin"+minsvname,sep='\t',index=False,\
        float_format='%5.8f')
	xmin2, ymin2, remG   = df2.nH, df2.Eadsavg, df2.Remarks
	df2.to_csv("GibbsFenemin"+minsvname,sep='\t',index=False,\
        float_format='%5.8f')
	#vibTScorre = dZPE - TScorr
	#ymin2 = ymin2 - TScorr +eU
	ymin2 = ymin2 +eU
	ax2.plot(x,[0.1 for i in range(len(x))],"--",color="k")
	ax2.plot(x,[0 for i in range(len(x))],"--",color="k")
	ax2.plot(x,[-0.1 for i in range(len(x))],"--",color="k")
	ax1.plot(x,y,'o')
	ax2.plot(x,y1,'o')
	ax1.plot(xmin,ymin,'--o')	
	ax2.plot(xmin2,ymin2,'--o')	
	#if input_file==Open1: # UNCOMMENT THE FOLLOWING LINES TO FIT THE CURVE
	#	from scipy.optimize import curve_fit as cfit 
	#	def model(x,a,b,c,d,e):
	#		y = c - a**(-((e*x)-b)/d)	
	#		return y
	#	def model2(x,a,b,c,d):
	#		y = a/(b+(d*np.exp(-(x-c)))) 
	#		return y
	#	xd = xmin2[4:]
	#	yd = ymin2[4:]
	#	popt,pcov = cfit(model,xd,yd)
	#	popt2,pcov2 = cfit(model2,xd,yd)
	#	x_extpl = np.linspace(5,60,1000)
	#	ax2.plot(x_extpl,model(x_extpl,*popt),alpha=0.3,linewidth=7) #FOR FITTING THE CURVE UNCOMMENT THIS LINE
	for i, txt in enumerate(rem):
		xytxt=(20,-20)
		txtsplt  = txt.split("#")
		txt = "".join(txtsplt)
		ax1.annotate(txt,(xmin[i],ymin[i]),fontsize=6,\
		textcoords="offset points",rotation=90,\
		xytext=xytxt,ha='center',va='bottom',\
		arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2'))
	for i, txtG in enumerate(remG):
		ax2.annotate(txtG,(xmin2[i],ymin2[i]),fontsize=6,\
		textcoords="offset points",rotation=90,\
		xytext=xytxt,ha='center',va='bottom',\
		arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2'))
	ax1.set_xlabel("Number of Hydrogen adsorbed")
	ax1.set_ylabel("Adsorption Energy (eV)")
	ax2.set_xlabel("Number of Hydrogen adsorbed")
	ax2.set_ylabel("Gibbs Free Energy per H (eV)")
	#ax2.set_xlim([x[0],x[-1]])
	fig.suptitle(system)
	fig.savefig(svname,format="png",dpi=500)
	#plt.show()
	return

def plotdGibbs(df,input_file,system):
	fname  = input_file.split(".")
	svname = "".join(fname[:-1])+pltDGib
	x = df['nH']
	y = df["diffGibbs"]
	plt.figure(figsize=(15,10))
	plt.plot(x,[0.1 for i in range(len(x))],"--",color="k")
	plt.plot(x,[0 for i in range(len(x))],"--",color="k")
	plt.plot(x,[-0.1 for i in range(len(x))],"--",color="k")
	df1 = df.sort_values(by=["nH","diffGibbs"])
	df2 = df1.drop_duplicates(subset=["nH"],\
        keep='first',ignore_index=True)
	xmin,ymin,rem = df2["nH"], df2["diffGibbs"],df2["Remarks(n,n-1)"]
	plt.plot(xmin,ymin,'--o')
	for i, txt in enumerate(rem):
		xytxt=(20,-20)
		txtsplt  = txt.split("#")
		txt = "".join(txtsplt)
		plt.annotate(txt,(xmin[i],ymin[i]),fontsize=6,\
		textcoords="offset points",rotation=90,\
		xytext=xytxt,ha='center',va='bottom',\
		arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2'))
	plt.xlabel("Adding nth Hydrogen")
	plt.ylabel("Differential Gibbs Free energy (eV)")
	plt.title(system+'\nThe labels represents n and n-1 position of H')
	plt.savefig(svname,format="png",dpi=500)
	return

def plotdAdsEne(df,input_file,system):
	fname  = input_file.split(".")
	svname = "".join(fname[:-1])+pltDAdsE
	x = df['nH']
	#vibTScorre = dZPE - TScorr
	#y = df["diffGibbs"] + TScorr - eU 
	y = df["diffGibbs"]  - eU 
	plt.figure(figsize=(15,10))
	plt.plot(x,y,'-o')
	df1 = df.sort_values(by=["nH","diffGibbs"])
	df2 = df1.drop_duplicates(subset=["nH"],keep='first',ignore_index=True)
	#vibTScorre = dZPE - TScorr
	#df2["diffGibbs"] = df2["diffGibbs"] + TScorr - eU
	df2.loc[:,"diffGibbs"] = df2.loc[:,"diffGibbs"]  - eU
	xmin,ymin,rem = df2["nH"], df2["diffGibbs"],df2["Remarks(n,n-1)"]
	for i, txt in enumerate(rem):
		xytxt=(20,-20)
		txtsplt  = txt.split("#")
		txt = "".join(txtsplt)
		plt.annotate(txt,(xmin[i],ymin[i]),fontsize=6,\
		textcoords="offset points",rotation=90,\
		xytext=xytxt,ha='center',va='bottom',\
		arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2'))
	plt.xlabel("Adding nth Hydrogen")
	plt.ylabel("Differential Adsorption Free energy (eV)")
	plt.title(system+'\nThe labels represents n and n-1 position of H')
	plt.savefig(svname,format="png",dpi=500)
	#plt.show()
	return


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#DEFINE CONSTANTS OF CALCULATION
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = 298 # K. This is the simulation temperature.
EH2        = -6.77026428   #IVDW-11 VALUE #eV
EslabAA7   = -781.49477915 #IVDW-11 VALUE #eV
EslabA7P   = -808.35898113 #IVDW-11 VALUE #eV 
eU = 0 #eV value for potential
EH2vib = 0.2732830  #IVDW-11 VALUE #eV THIS IS VIB-ENT ENERGY AT 298 K
entrEH2= 0.00135317888 # eV/K per H2 at 298 K JANAF Thermochemical tables 2ed. 
EH2tot = EH2 + EH2vib - (298)*entrEH2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#DEFINE FILES AND STORAGE TO BE USED
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Open1    = "AA7cvg.dat"
Open2 	 = "AA7+Pcvg.dat"
saveAvg  = "AvgEads.dat"
saveGibbs= "GibbsFene.dat"
savedGib = "diffGibbs.dat"
pltEavg  = "AvgEads.png"
pltDAdsE = "diffAdsEne.png"
pltDGib  = "diffGibbs.png" 
sysAA7   = "7 layer - Ni3P2 termination\nHydrogen adsorption"
sysAA7P  = "7 layer - Ni3P2+P termination\nHydrogen adsorption"
VibAA7   = "AA7VibEneSite.dat"
VibAA7P  = "AA7+PVibEneSite.dat"
col1 = "VibEne"
col2 = "Remarks"
dictAA7  = {}
dictAA7P = {}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# CALL FUNCTIONS BELOW
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dfAA7,dictAA7   = VibEneCollect(VibAA7,dictAA7)
dfAA7P,dictAA7P = VibEneCollect(VibAA7P,dictAA7P)

write_AvgEads(Open1,EslabAA7,sysAA7,VibAA7)
write_diffGibbs(Open1,EslabAA7,sysAA7,VibAA7)

write_AvgEads(Open2,EslabA7P,sysAA7P,VibAA7P)
write_diffGibbs(Open2,EslabA7P,sysAA7P,VibAA7P)

write_vib_entropy(VibAA7)
write_vib_entropy(VibAA7P)




