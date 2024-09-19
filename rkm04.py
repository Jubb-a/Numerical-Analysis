import math
import sympy as sym
from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tikzplotlib as tk

# import os

# os.chdir('C:\\Users\\ajubb\\Documents\\personal\\Math 501B\\project')
# df = pd.read_excel(r"C:\Users\ajubb\OneDrive - UCI Health\Documents\personal\Math 501A\De-identified Pt Fall Data for Analysis.xlsx")
df = pd.read_excel('De-identified Pt Fall Data for Analysis-copy.xlsx')
df_array = pd.DataFrame.to_numpy(df)

x,y,t= sym.symbols(('x','y','t'))
a=0.07 #damping coeff
b=5 #spring constant
f1=y
f2=(-b*x)-(a*y)
def RKMO4_IVP(f1, f2, x10, x20, h, t0, tn):
    tp = np.arange(t0+1.25,tn+h,h)
   # tprev = []
   # for i in range(len(tp)):
    #    tprev.append(tp[-i-1])
    n=len(tp)
    x1 = [x10]
    y2= [x20]
    for i in range(1,n):
        k11 = h*sym.N(f1.subs([(x, x1[i - 1]), (y, y2[i - 1]), (t,tp[i-1])]))
        k12 = h*sym.N(f2.subs([(x, x1[i - 1]), (y, y2[i - 1]),(t, tp[i-1])]))
        k21 = h*sym.N(f1.subs([(x, x1[i - 1] + 1/2*k11),(y, y2[i - 1] + 1/2*k12),(t, tp[i-1] + h/2)]))
        k22 = h*sym.N(f2.subs([(x, x1[i - 1] + 1/2*k11), (y, y2[i - 1] + 1/2*k12),(t, tp[i-1] + h/2)]))
        k31 = h*sym.N(f1.subs([(x, x1[i - 1] + 1/2*k21), (y, y2[i - 1] + 1/2*k22),(t, tp[i-1] + h)]))
        k32 = h*sym.N(f2.subs([(x, x1[i - 1] + 1/2*k21), (y, y2[i - 1] + 1/2*k22),(t, tp[i-1] + h)]))
        k41 = h*sym.N(f1.subs([(x, x1[i - 1] + k31), (y, y2[i - 1] + k32),(t, tp[i-1] + h)]))
        k42 = h*sym.N(f2.subs([(x, x1[i - 1] + k31), (y, y2[i - 1] + k32),(t, tp[i-1] + h)]))
        x1.append(x1[i - 1] + (k11 + 2*k21 + 2*k31 + k41)/6)
        y2.append(y2[i - 1] + (k12 + 2*k22 + 2*k32 + k42)/6)
    plt.plot(tp,x1,'c',df_array[:27,5],'ko-')
    plt.ylim(-0.8,0.8)
    plt.xlabel("Month")
    plt.ylabel('Patient Fall Rate (%)')
    #plt.legend(['RKMO4 Approximation', 'CNL Data'], loc='best')
    tex=tk.save(r'C:\Users\ajubb\OneDrive - UCI Health\Documents\personal\Math 501B\project\RKMO4 figure 1.tex')
    # tex=tk.get_tikz_code()
    return tp, x1, y2 #, tex
