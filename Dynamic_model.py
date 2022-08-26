# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
#import sympy as sp

#DH ---> Array de variables alpha, a, b, theta

class NE: 
    """
    Calculo de modelo dinamica mediante newton euler
    """
    def __init__(self,DH, s, m, I): #constructor
        self.DH = DH
        self.alpha= self.DH[:,0]
        self.a= self.DH[:,1]
        self.d= self.DH[:,2]
        self.theta= self.DH[:,3]
        ###########################
        self.p = np.transpose(np.array([[self.a],[self.d*(np.sin(self.alpha))],[self.d*(np.cos(self.alpha))]]))
        ###########################:,1
        self.DHshape = np.shape(self.DH)
        self.dof = self.DHshape[0]
        self.s = s
        self.m = m
        self.I = I

    def get_R(self):
        """
        esta funcion optiene la matrix de rotacion de cada articulacion 
        """
        theta, alpha, dof = self.theta,self.alpha, self.dof
        #self.d,self.a,
    
        R = [0]*(dof) #crear matriz de ceros de tamaño de los DOF
        for i in range(dof):
            ct = np.cos(theta[i])
            st = np.sin(theta[i])
            ca = np.cos(alpha[i])
            sa = np.sin(alpha[i])
            R[i]=np.array([[ct,-ca*st,sa*st],[st,ca*ct,-sa*ct],[0,sa,ca]])
        return R
    
    def get_Rt(self):
        dof = self.dof
        R2 = [0]*(dof)
        R = self.get_R()
        for i in range(dof):
            R2[i] = R[i].transpose()        
        return R2
    """
    forwards
    """
    def get_w(self,q_dot):
        dof = self.dof
        Ri = self.get_Rt()  #sacar la matriz de rotacion para casa Si
        w=[0]*(dof)  #crear una lista de ceros
        z0 = np.array([0, 0, 1])
        for i in range(dof):
            if i == 0:
                w_i = np.array([0, 0, 0])
                qq = q_dot[i]
            else:
                w_i = w[i-1]
                qq = q_dot[i]
            w[i] = Ri[i].dot(w_i+(z0*qq))
            #print(w[i])
        return w
    def get_w_dot(self,q_dot,q_2dot):
        dof = self.dof
        Ri = self.get_Rt()   #matrices de rotacion de cada Si
        w = self.get_w(q_dot)  #
        w_dot=[0]*(dof)  #creo lista 
        z0 = np.array([0, 0, 1])
        for i in range(dof):
            if i == 0:
                qqq= q_2dot[i]
                qq = q_dot[i]
                w_i = np.array([0, 0, 0])
                w_dot_i = np.array([0, 0, 0])
            else:
                qqq= q_2dot[i]
                qq = q_dot[i]
                w_i = w[i-1]
                w_dot_i = w_dot[i-1]
            w_dot[i]= Ri[i].dot(w_dot_i+(z0*qqq))+np.cross(w_i,(z0*qq))
            #print(w_dot[i])
        return w_dot
    
    def get_a(self,q_dot,q_2dot):
        dof = self.dof
        Ri = self.get_Rt()
        w = self.get_w(q_dot)
        w_dot = self.get_w_dot(q_dot,q_2dot)
        v_dot=[0]*(dof)
        for i in range(dof):
            if i==0:
                v_dot_i = np.array([0, 0, 0]) #velicidad de la base es cero
            else:
                v_dot_i = v_dot[i-1]

            w_i = w[i]
            w_dot_i = w_dot[i]
            pi = self.p[i]  
            v_dot[i]= np.cross(w_dot_i,pi)+ np.cross(w_i,np.cross(w_i,pi))+Ri[i].dot(v_dot_i)
            v_dot[i]=v_dot[i][0]
            #print(v_dot[i])
        return v_dot
    
    def get_acm(self,q_dot,q_2dot):
        dof = self.dof
        w = self.get_w(q_dot)
        w_dot = self.get_w_dot(q_dot,q_2dot)
        v_dot = self.get_a(q_dot,q_2dot)
        s = self.s
        acm=[0]*(dof)
        for i in range(dof):
            v_dot_i = v_dot[i]
            w_i = w[i]
            w_dot_i = w_dot[i] 
            si = s[i]  
            acm[i]= np.cross(w_dot_i,si)+ np.cross(w_i,np.cross(w_i,si))+v_dot_i
        return acm 
    """
    Backwards
    """
    def get_f(self,q_dot,q_2dot):
        R = self.get_R()  #la original, no la transpuesta
        a = self.get_acm(q_dot,q_2dot)
        m = self.m
        dof = self.dof
        f=[0]*(dof)
        for i in range(dof-1,-1,-1):
            if i == (dof-1): 
                fi= np.array([0, 0, 3]) #fuerza gravitacional de 0.3kg  
            else:
                fi = f[i+1]
            Ri = R[i]
            mi = m[i]
            ai = a[i]
            
            f[i]= Ri.dot(fi)+mi*ai
        return f
    def get_n(self,q_dot,q_2dot):
        R = self.get_R()  #la original no la transpuesta
        Rt = self.get_Rt()
        s = self.s
        dof = self.dof
        w = self.get_w(q_dot)
        w_dot = self.get_w_dot(q_dot,q_2dot)
        a = self.get_acm(q_dot,q_2dot)
        f = self.get_f(q_dot,q_2dot)
        m = self.m
        I = self.I
        p = self.p
        n=[0]*(dof)
        for i in range(dof-1,-1,-1):
            if i == (dof-1): 
                ni = np.array([0, 0, 0]) #sin torque en el efector final
                fi = np.array([0, 0, 3])
            else:
                ni = n[i+1] 
                fi = f[i+1]
            Ri = R[i]
            Rti = Rt[i]
            mi = m[i]
            ai = a[i]
            wi = w[i]
            w_dot_i = w_dot[i]
            si = s[i]
            pi = p[i]
            Ii = I[i]
            n[i] = Ri.dot(ni+np.cross(Rti.dot(pi[0]),fi))+np.cross((pi+si),(mi*ai))+Ii.dot(w_dot_i)+np.cross(wi,(Ii.dot(wi)))
            n[i]=n[i][0]
        return n
    
    def get_T(self,q_dot,q_2dot):
        Rt = self.get_Rt()
        dof = self.dof
        z0 = np.array([0, 0, 1])
        n = self.get_n(q_dot,q_2dot)
        T=[0]*(dof)
        for i in range(dof):
            T[i]=n[i].dot(Rt[i].dot(z0))
        return T
        
q1=1
q2=2
DH1 = np.array([[0,0,10,q1],[0,10,0,q2]])
s = np.array([[5,0,0],[5,0,0]])
m = [1,1]
I = [np.array([[1,0,0],[0,1,0],[0,1,0]]),np.array([[1,0,0],[0,1,0],[0,1,0]]) ]
CN = NE(DH1,s,m,I)
CN1 = NE(DH1,s,m,I)
"""
jj = CN.get_R()
ii = CN.get_Rt()
w = CN.get_w([1,1])
w_dot = CN.get_w_dot([1,1],[1,1])
v_dot = CN.get_a([1,1],[1,1])
acm = CN.get_acm([1,1],[1,1])
f = CN.get_f([1,1],[1,1])
n = CN.get_n([1,1],[1,1])
"""
T = CN.get_T([1,1],[1,1])
print(T)
