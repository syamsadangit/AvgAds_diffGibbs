import pandas as pd
import numpy as np

VibAA7  = "AA7VibEneSite2.dat"
VibAA7P = "AA7+PVibEneSite2.dat"

col1 = "VibEne"
col2 = "Remarks"

dictAA7  = {}
dictAA7P = {}

def VibEneCollect(input_file,dict):
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #READS INPUT FILE CONTAINIG VIB ENERGY AS DF AND DICT
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    col1 = "VibEne"
    col2 = "Remarks"
    with open(input_file) as file:
        next(file)
        for line in file:
            dat = [elem.rstrip() for elem in line.split()]
            key = dat.pop().split("#")[1] #removed last element from dat
            dict[key] = np.array([round(float(i), 3) for i in dat])
    rem = dict.keys()
    ene = dict.values()
    df = pd.DataFrame(list(zip(ene, rem)), columns=[col1, col2])
    print(" ")
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("Vibrational Energy Colleced from {}".format(input_file))
    print(df)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print(" ")
    return df,dict



df1, dictAA7  = VibEneCollect(VibAA7,dictAA7)
df2, dictAA7P = VibEneCollect(VibAA7P,dictAA7P)

#b = df1.loc[2,"VibEne"][2]
#a = type(b)
#print(b)
#print(a)
