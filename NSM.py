import numpy as np
import sympy as sym
import math
from numpy import linalg
import matplotlib.pyplot as plt

def f(x,y,yp):
    f = -y-yp
    return f

def fy(xp,z,zp):
    fy = z
    return fy

def fyp(xp,zp,zpp):
    fyp = zp
    return fyp

def NSM(alpha, beta, a, b, n, tol, maxiter):
    h=(beta-alpha)/n
    tk=(beta-alpha)/(b-a)
    count=1
    while count<=maxiter:
        x1=[alpha]
        x2=[tk]
        u1=0
        u2=1
        for i in range(1,n):
            x= a+(i-1)*h
            k11=h*x2[i-1]
            k12=h*f(x,x1[i-1],x2[i-1])
            k21=h*(x2[i-1]+(k12/2))
            k22=h*f((x+(h/2)),x1[i-1]+(k11/2),x2[i-1]+(k12/2))
            k31=h*(x2[i-1]+(k22/2))
            k32=h*f(x+h,x1[i-1]+(k21/2),x2[i-1]+(k22/2))
            k41=h*(x2[i-1]+k32)
            k42=h*f(x+h,x1[i-1]+k31,x2[i-1]+k32)
            x1.append(x1[i-1]+((k11+(2*k21)+(2*k31)+k41)/6))
            x2.append(x2[i-1]+((k12+(2*k22)+(2*k32)+k42)/6))
            kp11=h*u2
            kp12=h*(fy(x,x1[i-1],x2[i-1])*u1)
            kp21=h*(u2+(kp12/2))
            kp22=h*((fy(x+(h/2),x1[i-1],x2[i-1])*(u1+(kp11/2)))+((fyp(x+(h/2),x1[i-1],x2[i-1]))*(u2+(kp12/2))))
            kp31=h*(u2+(kp22/2))
            kp32=h*((fy(x+(h/2),x1[i-1],x2[i-1])*(u1+(kp21/2)))+((fyp(x+(h/2),x1[i-1],x2[i-1]))*(u2+(kp22/2))))
            kp41=h*(u2+kp32)
            kp42=h*((fy(x+h,x1[i-1],x2[i-1])*(u1+(kp11/2)))+((fyp(x+(h/2),x1[i-1],x2[i-1]))*(u2+(kp32/2))))
            u1= u1+((kp11+(2*kp21)+(2*kp31)+kp41)/6)
            u2= u2+((kp12+(2*kp22)+(2*kp32)+kp42)/6)
            r=abs(x1[i]-beta)
            if r<=tol:
                for i in range(len(x1+1)):
                    x=a+i*h
            tk=tk-((x1[-i]-beta)/u1)
            count+=1
            if count>= maxiter:
                print(count,'maxhit')
                break
    return x,x1,x2, plt.plot(x1,x2)

