import numpy as np
import sqlite3


class ReversibleRxn():
    def __init__(self, rxn):
        self.s=rxn.specieslist
        self.r=rxn.r_stoich
        self.p=rxn.p_stoich
        self.x=rxn.conc_list
        self.nuij=self.p - self.r
        self.kf=np.array(rxn.k_list)
        self.p0 = 1.0e+05 
        self.R = 8.3144598 
        self.gamma = np.sum(self.nuij, axis=0)
        self.length=len(self.s)

    def read_SQL(self, T):

        def chose_t_range(T):
            t_range=[]
            for species in self.s:
                v=cursor.execute('''SELECT THIGH
                from LOW WHERE species_name= ?''',(species,)).fetchall()  
                if v[0] > T:
                    t_range.append('high')
                else:
                    t_range.append('low')
            return t_range

        def get_coeffs(species_name, temp_range):
            if temp_range == 'low':
                v=cursor.execute('''SELECT COEFF_1,COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7 
                from LOW WHERE species_name= ?''',(species_name,)).fetchall()
            elif temp_range == 'high':
                v=cursor.execute('''SELECT COEFF_1,COEFF_2,COEFF_3,COEFF_4,COEFF_5,COEFF_6,COEFF_7 
                from HIGH WHERE species_name= ?''',(species_name,)).fetchall()
            coeffs=v[0]
            return coeffs 

        db = sqlite3.connect('NASA.sqlite') 
        cursor = db.cursor()
        coefs=[]
        t_range=chose_t_range(T)
        s_t=zip(self.s, t_range)
        for species, tmp in s_t:
            coef=get_coeffs(species,tmp)
            coefs.append(coef)
        self.nasa=np.array(coefs)

  

    def Cp_over_R(self, T):
        a = self.nasa
        Cp_R = (a[:,0] + a[:,1] * T + a[:,2] * T**2.0 
                + a[:,3] * T**3.0 + a[:,4] * T**4.0)
        return Cp_R

    def H_over_RT(self, T):
        a = self.nasa
        H_RT = (a[:,0] + a[:,1] * T / 2.0 + a[:,2] * T**2.0 / 3.0 
                + a[:,3] * T**3.0 / 4.0 + a[:,4] * T**4.0 / 5.0 
                + a[:,5] / T)
        return H_RT
               
    def S_over_R(self, T):
        a = self.nasa
        S_R = (a[:,0] * np.log(T) + a[:,1] * T + a[:,2] * T**2.0 / 2.0 
               + a[:,3] * T**3.0 / 3.0 + a[:,4] * T**4.0 / 4.0 + a[:,6])
        return S_R

    def backward_coeffs(self, T):

        # Change in enthalpy and entropy for each reaction
        delta_H_over_RT = np.dot(self.nuij.T, self.H_over_RT(T))
        delta_S_over_R = np.dot(self.nuij.T, self.S_over_R(T))

        # Negative of change in Gibbs free energy for each reaction 
        delta_G_over_RT = delta_S_over_R - delta_H_over_RT

        # Prefactor in Ke
        fact = self.p0 / self.R / T

        # Ke
        ke = fact**self.gamma * np.exp(delta_G_over_RT)

        self.kb = self.kf / ke

    def prograss_rate(self,T):

        self.read_SQL(T)
        self.backward_coeffs(T)

        omega=self.kf*np.product(self.x**self.r.T)-self.kb*np.product(self.x**self.p.T)
        return omega

    def rxn_rate(self,T):
        omega=self.prograss_rate(T)
        return np.sum(omega * self.nuij, axis=1)


