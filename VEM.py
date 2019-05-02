import numpy as np
from scipy.special import gamma,psi,polygamma

class VEM:

    def __init__(self, h, var_h, fe, L_g, a):

        self.eps = np.finfo(float).eps

        self.h = h # Observed variable : RIR
        self.var_h = var_h # variance du bruit de mesure
        self.fe = fe # sample rate
        self.L_h_true = h.shape[0] # length of h
        self.L_g = L_g # length of microphone+source IR
        self.L_h = self.L_h_true - self.L_g + 1
        self.a = a # Absortion coefficient


        self.g = np.zeros(self.L_g) + 1
        self.alpha = np.zeros(self.L_h) + 1
        self.beta = np.zeros(self.L_h) + 1

        self.lamb = 1
        self.b = np.zeros(self.L_h_true) + 1
        self.h_hat = np.zeros(self.L_h_true) + 1

        self.u = np.arange(self.L_h_true)/self.fe
        print(self.u.shape)
        self.v = np.arange(self.L_h)/self.fe + self.eps

    def vfe(self):

        L_h_term = self.L_h_true/2*np.log(2*np.pi*self.var_h)
        var_h_term = np.exp(-2*self.a*self.u)/(2*self.var_h)
        quadr_term = (np.exp(self.a*self.u)*np.squeeze(self.h) - np.convolve(self.g,self.alpha/self.beta,mode='full'))**2
        var_term = np.convolve(self.g**2,self.alpha/(self.beta)**2, mode='full')
        likelihood = - L_h_term - np.sum(var_h_term*(quadr_term + var_term))

        temp1 = self.lamb*self.v**2*np.log(self.v) + np.log(gamma(self.alpha)/gamma(self.lamb*self.v**2))
        temp2 = (self.lamb*self.v**2 - self.alpha)*(psi(self.alpha) - np.log(self.beta)) - self.alpha*(np.log(self.beta + self.v/self.beta - 1))
        pi_term = np.sum(temp1 + temp2)

        return likelihood + pi_term

    def update_alpha(self):

        for v in range(self.L_h):
            temp1 = np.sum()
