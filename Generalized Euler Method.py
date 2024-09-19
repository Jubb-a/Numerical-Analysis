#Math 501B Lecture 1 Example 2
import math
import numpy as np
import matplotlib.pyplot as plt

def myEM(f,t0,tn,h,x0): #takes f as lambda function
    h=h
    t= np.arange(t0,tn+h,h)
    w=[x0]
    for i in range(len(t)-1):
        wn = w[i] + h*f(w[i], t[i])
        w.append(wn)
    return t,w, plt.plot(t,w)
    
