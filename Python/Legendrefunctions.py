# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:18:35 2019

@author: Amin Khiabani
 aminkhiabani@outlook.com
"""
import math

def factorial(n) :    
    if(n == 0):
       return  1
    else:
      return n*factorial(n - 1)
         

def legendf(x,n) :    
    if(n==0):
        return 1.0
    elif(n==1):
        return x
    else:
        return (-(n-1)*legendf(x,n-2) + (2*n-1)*x*legendf(x,n-1))/n
    
def legendrefrange(x,n):
    
    LFS = []
        
    LFS.append(1.0)
    LFS.append(x)
    
    if( n == 0):
        return LFS[0]
    elif ( n == 1) :
        return LFS
    else:    
        for I in range(2,n):         
            LFS.append((-(I-1)*LFS[I-2] + (2*I-1)*x*LFS[I-1])/I)
        return LFS
    
def Legendre(x,n):    
    
    LFS0 = 1.0
    LFS1 = x
    Legen = 0.0
    if(n==0) :
        return LFS0                  
    elif(n==1) :
        return LFS1
    elif(n>1):        
        for I in range(2,n+1):
            Legen = (-(I-1)*LFS0 + (2*I-1)*x*LFS1)/(I)
            LFS0 = LFS1
            LFS1 = Legen
    return Legen
        
    
def AlegendreRange(x,n):
     
    s = math.sqrt(1.0-x*x)
     
    Pnm =[]
     
    if(n==0):
         Pnm.append(1.0)
    elif(n==1) :
         Pnm.append(x)
         Pnm.append(s)
    elif(n==2) :           
         Pnm.append((3.0*x*x-1.0)/2.0)
         Pnm.append(3.0*s*x)
         Pnm.append(3.0*(1.0-x*x))
    elif(n==3) :
        Pnm.append(x*(5.0*x*x-3.0)/2.0)
        Pnm.append(3.0*(5.0*x*x-1.0)*s/2.0)
        Pnm.append(15.0*x*(1.0-x*x))
        Pnm.append(15.0*s*s*s)
    elif(n==4) :
        Pnm.append((35.0*x*x*x*x-30.0*x*x+3.)/8.0)
        Pnm.append(5.0*(7.0*x*x-3.0)*x*s/2.0)
        Pnm.append(15.0*(7.0*x*x-1.0)*(1.0-x*x)/2.0)
        Pnm.append(105.0*s*s*s*x)
        Pnm.append(105.0*s*s*s*s)
    elif(n==5) :
        Pnm.append(x*(63.0*x*x*x*x-70.0*x*x+15.0)/8.0)
        Pnm.append(15.0*s*(21.0*x*x*x*x-14.0*x*x+1.)/8.0)
        Pnm.append(105.0*x*(3.0*x*x-1.0)*(1.0-x*x)/2.0)
        Pnm.append(105.0*s*s*s*(9.*x*x-1.)/2.0)
        Pnm.append(945.0*s*s*s*s*x)
        Pnm.append(945.0*s*s*s*s*s)
    elif(n==6) :
        Pnm.append((x*x*(x*x*(231.0*x*x-315.0)+105.0)-5.0)/16.0)
        Pnm.append(21.*x*(x*x*(33.0*x*x-30.0)+5.0)*s/8.0)
        Pnm.append(105.0*s*s*(x*x*(33.0*x*x-18.0)+1.)/8.0)
        Pnm.append(315.0*(11.0*x*x-3.0)*x*s*s*s/2.0)
        Pnm.append(945.0*s*s*s*s*(11.0*x*x-1.0)/2.0)
        Pnm.append(10395.0*x*s**5)
        Pnm.append(10395.0*s**6)
    elif(n==7) :
        Pnm.append(x*(x*x*(429.0*x**4-693.0*x*x+315.0)-35.0) /16.0)
        Pnm.append(7.0*s*(x*x*(429.0*x**4-495.0*x*x+135.0)-5.0)/16.0)
        Pnm.append(63.0*x*s*s*(x*x*(143.0*x*x-110.0)+15.0)/8.0)
        Pnm.append(315.0*s*s*s*(x*x*(143.0*x*x-66.0)+3.)/8.0)
        Pnm.append(3465.0*x*s*s*s*s*(13.0*x*x-3.0)/2.0)
        Pnm.append(10395.0*(s**5)*(13.0*x*x-1.0)/2.0)
        Pnm.append(135135.0*x*s**6)
        Pnm.append(135135.0*s**7)
    elif(n>=8) :
        ALFS = []
        ALFS.append((x*x*(x*x*(x*x*(6435.0*x*x-12012.0)+6930.0)-1260.0)+35.0)/128.0)
        ALFS.append(9.0*x*s*(x*x*(x*x*(715.0*x*x-1001.0)+385.0)-35.0)/16.0)
        ALFS.append(315.0*s*s*(x*x*(x*x*(143.0*x*x-143.0)+33.0)-1.0)/16.0)
        ALFS.append(3465.0*x*s*s*s*(x*x*(39.0*x*x-26.0)+3.0)/8.0)
        ALFS.append(10395.0*s*s*s*s*(65.0*x**4-26.0*x*x+1.0)/8.0)
        ALFS.append(135135.*(x*s**5)*(5.0*x*x-1.0)/2.0)
        ALFS.append(135135.0*(s**6)*(15.0*x*x-1.)/2.0)
        ALFS.append(2027025.0*x*s*s*s*s*s*s*s)
        ALFS.append(2027025.0*s*s*s*s*s*s*s*s)
        
        if(n == 8) :
            Pnm = ALFS
        elif(n>8) :  
            BLFS = []
            BLFS.append(x*(x*x*(429.0*x**4-693.0*x*x+315.)-35.0) /16.0)
            BLFS.append(7.0*s*(x*x*(429.0*x**4-495.0*x*x+135.0)-5.0)/16.0)
            BLFS.append(63.0*x*s*s*(x*x*(143.0*x*x-110.0)+15.0)/8.0)
            BLFS.append(315.0*s*s*s*(x*x*(143.0*x*x-66.0)+3.0)/8.0)
            BLFS.append(3465.0*x*s*s*s*s*(13.0*x*x-3.0)/2.)
            BLFS.append(10395.0*(s**5)*(13.0*x*x-1.0)/2.)
            BLFS.append(135135.0*x*s**6)
            BLFS.append(135135.0*s**7)
            
            l = 8
            while(l<n):                
                Pnm = []            
                Pnm.append(Legendre(x,l+1))
                k = 1
                while(k<=l+1):                                                        
                    if(k<=l-1) :
                        Pnm.append(((2*l+1)*x*ALFS[k]-(l+k)*BLFS[k])/(l-k+1))
                    elif(k>l-1) :
                        Pnm.append(-((l-k+2)*x*Pnm[k-1]-(l+k)*ALFS[k-1])/s)                                          
                    k = k + 1                                  
                BLFS = ALFS               
                ALFS = Pnm                                      
                l = l + 1
    return Pnm 

def AlegendreRangeTest(x,n):
     
    s = math.sqrt(1.0-x*x)
     
    Pnm =[]
     
    if(n==0):
         Pnm.append(1.0)
    elif(n==1) :
         Pnm.append(x)
         Pnm.append(s)
    elif(n>=2) :
        ALFS = []           
        ALFS.append((3.0*x*x-1.0)/2.0)
        ALFS.append(3.0*s*x)
        ALFS.append(3.0*(1.0-x*x))      
        if(n == 2) :
            Pnm = ALFS
        elif(n>2) :  
            BLFS = []
            BLFS.append(x)
            BLFS.append(s)
            
            l = 2
            while(l<n):                
                Pnm = []            
                Pnm.append(Legendre(x,l+1))
                k = 1
                while(k<=l+1):                                                        
                    if(k<=l-1) :                       
                        Pnm.append(((2*l+1)*x*ALFS[k]-(l+k)*BLFS[k])/(l-k+1))
                    elif(k>l-1) :                       
                        Pnm.append(-((l-k+2)*x*Pnm[k-1]-(l+k)*ALFS[k-1])/s)                                          
                    k = k + 1                                                  
                BLFS = ALFS               
                ALFS = Pnm                                      
                l = l + 1                 
    return Pnm                         
                       
                   
def FulNormALegenSerie(x,n):
      
      s = math.sqrt(1.0-x*x)
      Pnm = []    
      BLFS = []
      ALFS = []     

      BLFS.append(x*math.sqrt(3.0)) 
      BLFS.append(s*math.sqrt(3.0))

      ALFS.append(math.sqrt(5.0)*Legendre(x,2)) 
      ALFS.append(3.0*s*x*math.sqrt(10.0/6.0))
      ALFS.append(3.0*(1.0-x*x)*math.sqrt(10.0/24.0))

      if(n==0) :
          Pnm.append(1.0) 
      elif(n == 1) :            
            Pnm = BLFS
      elif(n==2):               
            Pnm= ALFS            
      elif(n>2):       
            l = 3
            while(l<= n):                  
                  Pnm = []
                  Pnm.append(math.sqrt(2.0*(l)+1.0)*Legendre(x,l)) 
                  k = 1
                  while(k <= l):
                        if(k < l-1) :
                              Anm = math.sqrt(((2*l-1)*(2*l+1))/((l-k)*(l+k)))

                              Bnm = -math.sqrt((((l-1)*(l-1)-k*k))/(4*(l-1)*(l-1)-1))

                              Pnm.append(Anm*(x*ALFS[k] + Bnm*BLFS[k]))

                        elif(k >= l-1) :
                              if(k == l-1):
                                    Pnm.append(math.sqrt(2.0*(k)+3.0)*x*ALFS[l-1])
                              elif(k == l) :                                    
                                    Pnm.append(s*math.sqrt(1.0+1.0/(2.0*k))*ALFS[l-1])                            
                        
                        k = k + 1
                  
                  BLFS = ALFS
                  ALFS = Pnm
                  l = l + 1                 
      
      return Pnm
  
def FulNormALegenRange(x,n):
      
      s = math.sqrt(1.0-x*x)
      Pnm = []     
      
      if(n <= 8) :
            Pnm = AlegendreRange(x,n)
            Pnm[0]= math.sqrt(2.0*(n)+1.0)*Legendre(x,n) 
            for k in range(1,n+1):
                ld = math.sqrt(2.0*((2*n+1)*factorial(n-k))/(factorial(n+k)))
                Pnm[k] = ld*Pnm[k]
            
      elif(n > 8) :            
            BLFS = AlegendreRange(x,7)
            for I in range(1,8):
                  ld = math.sqrt(30.0*(factorial(7-I))/(factorial(7+I)))
                  BLFS[I] = ld*BLFS[I]           
            
            ALFS= AlegendreRange(x,8)
            for I in range(1,9):
                  ld = math.sqrt(34.0*(factorial(8-I))/(factorial(8+I)))
                  ALFS[I] = ld*ALFS[I]
            
            l = 9
            while(l<= n):                  
                  Pnm = []
                  Pnm.append(math.sqrt(2.0*(l)+1.0)*Legendre(x,l)) 
                  k = 1
                  while(k <= l):
                        if(k < l-1) :
                              Anm = math.sqrt(((2*l-1)*(2*l+1))/((l-k)*(l+k)))

                              Bnm = -math.sqrt((((l-1)*(l-1)-k*k))/(4*(l-1)*(l-1)-1))

                              Pnm.append(Anm*(x*ALFS[k] + Bnm*BLFS[k]))

                        elif(k >= l-1) :
                              if(k == l-1):
                                    Pnm.append(math.sqrt(2.0*(k)+3.0)*x*ALFS[l-1])
                              elif(k == l) :
                                    #ld = float(k)
                                    Pnm.append(s*math.sqrt(1.0+1.0/(2.0*k))*ALFS[l-1])                            
                        
                        k = k + 1
                  
                  BLFS = ALFS
                  ALFS = Pnm
                  l = l + 1                   
      
      return Pnm
    
def gravitation(phi,Lambda,fileName):
    
    x = math.sin(phi*math.pi/180.0)
    Rlam = Lambda*math.pi/180.0 
        
    
    Sum1 = 0.0
    Tsum = 0.0
    
    try:    
        file = open(fileName,'r')
        for line in file:
            if(line.find('earth_gravity_constant')>=0):        
                l = line.rsplit()
                if(l[1].find('D') != -1):                    
                    l[1] = l[1].replace('D','e')                    
                G = float(l[1])
                print('earth_gravity_constant= {:20.4f}'.format(G))
            elif(line.find('radius')>=0):
                l = line.rsplit()
                if(l[1].find('D') != -1):                    
                    l[1] = l[1].replace('D','e')      
                R = float(l[1])
                print('radius= {:10.2f}'.format(R))
                r, Rn = ellipsoidal(R,1,phi) 
            elif(line.find('max_degree')>=0):
                l = line.rsplit()
                N = int(l[1])
                print('max_degree= {}'.format(N))                 
            elif(line.find('gfc')>=0) :
                l = line.rsplit()
                L = int(l[1])
                M = int(l[2])
                C = float(l[3])
                S = float(l[4])                
                if(M==0) :
                    Pn = FulNormALegenSerie(x,L)
                    Sum1 = Pn[M]*C
                else:
                    Sum1 = Sum1 + Pn[M]*(C*math.cos(M*Rlam)+ \
                            S*math.sin(M*Rlam))                    
                if(M==L):
                 
                     Tsum = Tsum + Sum1*(L-1)*Rn**L
                 
                  #   print("{0} , {1} , {2} , {3}".format(L,M,C,S))
    except IOError :
        print('Error in reading file')  
        file.close()
        
    file.close()
    return Tsum*G/(r*r)
    
def normalGravity(latitute): 

    e2 = 0.00669438002290
    k = 0.001931851353
    
    DEGRAD = math.pi/180.
    e4 = e2*e2
    e6 = e4*e2
    
    Rlat = latitute*DEGRAD
    x = math.sin(Rlat)
    Term= 0.0
    
    a2n = []
    a2n.append( 0.50*e2 + k)    
    a2n.append(3.0*e4/8.0 + 0.50*e2*k)    
    a2n.append(5.0*e6/16.0 + 3.0*e4*k/8.0)
    a2n.append(35.0*e4*e4/128.0 + 5.0*e6*k/16.0)
    i = 1
    for l in a2n:
            Term = Term + l*x**(2*i) 
            i = i + 1

    return 9.78032677150*(1.0+Term)

def adaptedGravity(Latitute,Height):
      
       ga = 9.78032677150     
       DEGRAD = math.pi/180.0

       x = math.sin(Latitute*DEGRAD)
       x2 = x*x

       return  ga*(1+x2*(0.00527904140 + x2*(0.00002327180+ \
            x2*(0.0000001262e0+0.0000000007e0*x2))))-Height*(0.03087798e-4 - \
            x2*( 0.0000439e-4 + 0.00000020e-4*x2))- Height*Height*(-0.00007265e-8 \
            + 0.00000021e-8*x2)  
        
def ellipsoidal(a,method,latitute):   

      PI = 3.1415926535897932384626430
      DEGRAD = PI/180.0
      e2 = 0.006694380022900
      ep2 = 0.006739496775480
      
      Rlat = latitute*DEGRAD
      x = math.sin(Rlat)
      x2 = x*x
      
      ellipsoid = []

      if(method == 0):
            ellipsoid.append(a/math.sqrt(1+ep2*x2))
      elif(method == 1) :
            w2 = 1.0-e2*x2
            ellipsoid.append(a*math.sqrt(1.0 + e2*(e2-2.0)*x2)/math.sqrt(w2))      
            
      ellipsoid.append(a/ellipsoid[0])
      
      return ellipsoid

def gravityCoef(method):
    
    import numpy as np
    
    a = 6378137.0  
    Omega = 0.7292115e-4
    GM = 0.39860050000e+15    
    DEGRAD = np.pi/180.0
    Omega2 = Omega*Omega    
    
    ga = [9.78032677150,normalGravity(30.0), \
          9.8061992030,normalGravity(60.0),9.8321863685]
   
    Phi = [0.0, 30.0, 45.0, 60.0 ,90.0]
    
    A = []
    B = []
    
    for i in range(len(Phi)):
        Rphi = Phi[i]*DEGRAD
        x = math.sin(Rphi)
        s = math.cos(Rphi)    
        s2 = s*s
        r,Rn = ellipsoidal(a,method,Phi[i])     
        B.append([ga[i]+ r*Omega2*s2])
        row = np.zeros(5,dtype=np.float64)
        for j in range(len(ga)):            
            row[j] = GM*((2*j+1)*Legendre(x,2*j)*Rn**(2*j+1))/(r*a)            
        A.append(row)
     
  #  B = np.matrix(B,dtype='f8')
  #  A = np.matrix(A,dtype='f8')  
  #  return (A**(-1))*B    
    A = np.array(A,dtype=np.float64)
    B = np.array(B,dtype=np.float64)
   
    return np.linalg.inv(A).dot(B)


def calcGravity(latitude,method):
    
    a = 6378137.0   
    Omega = 0.7292115e-4
    GM = 0.39860050000e+15
    
    DEGRAD = math.pi/180.0
    Omega2 = Omega*Omega
    
    Rphi = latitude*DEGRAD
    x = math.sin(Rphi)
    s = math.cos(Rphi)
    s2 = s*s
    
    r,Rn = ellipsoidal(a,method,latitude)  
       
    C = gravityCoef(method)
    
    Sigma = 0.0
    
    for j in range(len(C)):
        Sigma = Sigma + GM*((2*j+1)*C[j]*Legendre(x,2*j)*Rn**(2*j+1))
        
    return Sigma/(r*a) - r*Omega2*s2 


    
