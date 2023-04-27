import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random as rd
import re 
from sys import exit
import pickle

pd.options.mode.chained_assignment = None

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#DEFINE FUNCTIONS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def iseven(num):
	if (num % 2 == 0):
		condition = True
	else:
		condition = False
	return condition

def FreeEntVib(T,Evib):
	kb = 8.617333262e-05 #eV/K 
	B = 1/(kb*T)
	x = B*Evib
	exp = np.exp(2*x)
	termvib     = Evib + (2*Evib/(exp - 1))
	termentropy = (1/B)*( (2*x/(exp - 1)) - np.log(1-(1/exp)) )   
	F = termvib - termentropy
	return F,termvib,termentropy

def write_vib_entropy(Vibenefile):
	#T = 298  #K
	head_file = "Entropy_correction"
	dictinp = {}
	dfvib,dictvib = VibEneCollect(Vibenefile,dictinp)
	col = dfvib.columns
	Evib = dfvib[col[0]].copy() #SettingWithCopyWarning
	Free_energy_entropy,termvib,termentropy = FreeEntVib(T,Evib)
	#print(Free_energy_entropy)
	#print(Evib)
	dfvib.insert(1,'Entrene',termentropy,True)
	dfvib.insert(2,'tot_corr',Free_energy_entropy,True)
	#print(dfvib)
	fname  = Vibenefile.split(".")
	svname = head_file+"".join(fname[:-1])+".dat"
	dfvib.to_csv(svname,sep='\t',index=False,\
	float_format='%5.8f')
	return
	 
			
		

def EneCorrectVib(Vibenefile,dfEads,Eslab):
	#This function takes the energy of the 
	#system after optimisation, Energy of 
	#desorbed slab and add the corresponding 
	#site-vibrational energy.
	print("The vib ene data is taken from"+Vibenefile)
	PHlst = ["PH2","PH3"] # update this to more num if needed
	#T = 298 #K The entropy values are for STP
	dictinp = {}
	dfvib,dictvib = VibEneCollect(Vibenefile,dictinp) #CALL FUNCTION
	remkeys = dictvib.keys()
	#print(dictvib)
	for ky in remkeys:
		Evib = dictvib[ky]
		dictvib[ky],V,E =  FreeEntVib(T,Evib) #CALL FUNCTION
	#print(dictvib)
	#col = dfvib.columns
	#Evib = dfvib[col[0]]
	#Free_energy_entropy,termvib,termentropy = FreeEntVib(T,Evib)

	if Eslab == EslabA7P:
		Eslab = Eslab+dictvib["Ni3P2+4P"]
	elif Eslab == EslabAA7:
		Eslab = Eslab
	else:
		print("Check the slab energy and vibrational energies") 
	dfeads = dfEads.copy()
	#print("Before Corection ----> ",dfeads)
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
	#print("After Corection ----> ",dfeads)
	return dfeads,Eslab


def VibEneCollect(input_file,dict):
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
	#Calculates Avg Ads Energy
	#avgEads = (E_slab+H - E_slab - (1/2 * n * EH2) )
	#df1 is Avg ads energy (no correction)
	#df2 is Gibbs free energy per H (Vib-TS corrected)
	df1 = df.copy()
	y = (df1.Eads - Eslab - (df1.nH*EH2*0.5))/df1.nH
	dfeads = df.copy()
	df2,Eslab = EneCorrectVib(Vibenefile,dfeads,Eslab) #CALL FUNCTION
	y2 = (df2.Eads - Eslab - (df2.nH*EH2tot*0.5))/df2.nH
	df1.insert(2,'Eadsavg',y,True)
	df2.insert(2,'Eadsavg',y2,True)
	return df1,df2 

def diffGibbs(df1,df2,neibr):
	# Calculates the differential Gibbs free energy
	#DelGn=En - En-1 - 0.5*EH2 + 0.24 + e*U (U=0 here)
	#df1 - without ref ene val. df2 - with ref ene val 
	#use df1,df2 = cleandf(df) to obtain df1 and df2
	dropind = len(df2.index) - 1
	dfn     = df1.copy()
	dfn_1   = df2.copy()
	arn     = dfn.values.tolist()   
	arn_1   = dfn_1.values.tolist() 
	#vibTScorre = dZPE - TScorr
	#terms   = (-0.5*EH2tot) - TScorr + eU
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
	print(df3)
	if (rem_avail):
		drpind = df3.loc[ ~df3["Remarks"].isin(rem) ].index	
		df3.drop(drpind,inplace=True)
	print(df3)
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
	global rem
	fname  = input_file.split(".")
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
	xmin2,ymin2   = df2.nH, df2.Eadsavg
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
	if input_file==Open1:
		from scipy.optimize import curve_fit as cfit 
		def model(x,a,b,c,d,e):
			y = c - a**(-((e*x)-b)/d)	
			return y
		def model2(x,a,b,c,d):
			y = a/(b+(d*np.exp(-(x-c)))) 
			return y
		xd = xmin2[4:]
		yd = ymin2[4:]
		popt,pcov = cfit(model,xd,yd)
		popt2,pcov2 = cfit(model2,xd,yd)
		#p  = np.polyfit(xd,yd,7)
		#print("Parameters from fit: ",p)
		#fit_fn = np.poly1d(p)
		#print(xd)
		x_extpl = np.linspace(5,60,1000)
		#y_extpl = fit_fn(x_extpl)
		#ax2.plot(x_extpl,y_extpl,alpha=0.1,linewidth=7)
		ax2.plot(x_extpl,model(x_extpl,*popt),alpha=0.3,linewidth=7)
		#ax2.plot(x_extpl,model2(x_extpl,*popt2),alpha=0.2,linewidth=10)
	for i, txt in enumerate(rem):
		xytxt=(20,-20)
		txtsplt  = txt.split("#")
		txt = "".join(txtsplt)
		ax1.annotate(txt,(xmin[i],ymin[i]),fontsize=6,\
		textcoords="offset points",rotation=90,\
		xytext=xytxt,ha='center',va='bottom',\
		arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2'))
		ax2.annotate(txt,(xmin2[i],ymin2[i]),fontsize=6,\
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
	#plt.plot(x,y,'o')
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
	#plt.ylim([-1.2,0.7])
	plt.ylabel("Differential Gibbs Free energy (eV)")
	plt.title(system+'\nThe labels represents n and n-1 position of H')
	plt.savefig(svname,format="png",dpi=500)
	#plt.show()
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
EH2        = -6.77025412   #-6.77024788 #eV
EslabAA7   = -741.72008203 #-741.72017502  #eV
EslabA7P   = -767.48135606  #eV 
#vibTScorre = 0.24 #eV
#TScorr	   = -0.20 #eV T = 293.15K
#dZPE       = 0.04 #eV
eU = 0 #eV value for potential
EH2vib = 0.536823477 #eV THIS IS VIB-ENT ENERGY AT 298 K
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




