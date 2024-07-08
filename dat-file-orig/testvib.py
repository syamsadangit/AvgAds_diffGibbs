import pandas as pd

VibAA7  = "AA7VibEneSite.dat"
VibAA7P = "AA7+PVibEneSite.dat"

col1 = "VibEne"
col2 = "Remarks"

dictAA7  = {}
dictAA7P = {}

def VibEneCollect(input_file,col1,col2,dict):
	df = pd.read_csv(input_file,skiprows=1,sep=" ",\
		names=(col1,col2))
	for i in range(len(df[col2])):
		key  = df[col2][i].split("#")[1]
		dict["{}".format(key)] = df[col1][i]
	print("Vibrational Energy Colleced from {}".format(input_file))
	print(dict) 

VibEneCollect(VibAA7,col1,col2,dictAA7)
VibEneCollect(VibAA7P,col1,col2,dictAA7P)

