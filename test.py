import numpy as np

def FreeEntVib(T,Evib):
        kb = 8.617333262e-05 #eV/K
        B = 1/(kb*T)
        x = B*Evib
        exp = np.exp(2*x)
        termvib     = Evib + (2*Evib/(exp - 1))
        termentropy = (1/B)*( (2*x/(exp - 1)) - np.log(1-(1/exp)) )
        F = termvib - termentropy
        return F,termvib,termentropy

e = FreeEntVib(273.15,0.536823477)
print(e)
