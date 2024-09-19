#Math 501B Homework 1 Section 8.4 Problem 1CP
#Write and test a subprogram or procedure for the fifth-order Adams-Bashforth-Moulton
#method coupled with a fourth-order Runge-Kutta method (see Computer Problems 8.3.7-9,
#p. 548). Keep only the five most recent values of (t_i,x_i). Print statements should be
#included in the routine.
import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sym
 
x,t=sym.symbols(('x','t'))
def RKMO5andAdamBashforth(f,t0,tm,h,x0):
    x,t=sym.symbols(('x','t'))
    f=f
    h=h
    w=[x0]
    p=np.arange(t0,tm+h,h)
    #RKMO5 loop for w0,w1,w2,w3,and w4
    for i in range(5):
        f1=h*f.subs([(x,w[i]),(t,p[i])])
        f2=h*f.subs([(x,(w[i]+((1/2)*f1))),(t,(p[i]+ (h/2)))])
        f3=h*f.subs([(x,(w[i]+((1/2)*f2))),(t,(p[i]+ (h/2)))])
        f4=h*f.subs([(x,(w[i]+f3)),(t,(p[i]+ (h/2)))])
        f5=h*f.subs([(x,(w[i]+((1/27)*((7*f1)+(10*f2)+f4)),(t,(p[i]+ (h/2)))))])
        f6=h*f.subs([(x,(w[i]+((28/625)*f1)-((1/5)*f2)+((546/625)*f3)+((54/625)*f4)-((378/625)*f5))),(t,(p[i]+ (2*h/3)))])
        wn=w[i]+((1/24)*f1)+((5/48)*f4)+((27/56)*f5)+((125/336)*f6)
        w.append(wn)
    #Adam Bashforth loop
    for i in range(4,5):
        wn=w[i]+((h/720)*((1901)*f.subs([(x,w[i]),(t,p[i])])-(2274*f.subs([(x,w[i-1]),(t,p[i-1])]))+(2616*f.subs([(x,w[i-2]),(t,p[i-2])]))-(1274*f.subs([(x,w[i-3]),(t,p[i-3])]))+(251*f.subs([(x,w[i-4]),(t,p[i-4])]))))
        w.append(wn)
    return p,w, plt.plot(p,w)

