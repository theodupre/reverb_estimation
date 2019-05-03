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


        self.g = np.zeros(self.L_g) + 0.005
        # self.g[0] = 1

        self.alpha = np.ones(self.L_h)*0.001
        self.beta = np.ones(self.L_h)*0.5

        self.lamb = 0.5
        self.b = np.ones(self.L_h_true)*0.5
        self.h_hat = np.ones(self.L_h_true)*0.5

        self.u = np.arange(1,self.L_h_true+1)
        self.v = np.arange(1,self.L_h+1)

        self.e_u = np.exp(-self.a*self.u)
        self.e_2u = np.exp(-2*self.a*self.u)

    def vfe(self):

        L_h_term = self.L_h_true/2*np.log(2*np.pi*self.var_h)
        var_h_term = np.exp(-2*self.a*self.u)/(2*self.var_h)
        quadr_term = np.square(np.squeeze(self.h)/self.e_u - np.convolve(self.g,self.alpha/self.beta,mode='full'))
        var_term = np.convolve(self.g**2,self.alpha/np.square(self.beta), mode='full')
        likelihood = - L_h_term - np.sum(var_h_term*(quadr_term + var_term))
        # print(gamma(self.alpha)) # /(gamma(self.lamb*self.v**2) + self.eps)
        temp1 = self.lamb*self.v**2*np.log(self.v + self.eps) + np.log(gamma(self.alpha)/(gamma(self.lamb*self.v**2) + self.eps) + self.eps )
        temp2 = (self.lamb*self.v**2 - self.alpha)*(psi(self.alpha) - np.log(self.beta)) - self.alpha*(np.log(self.beta + self.v/self.beta - 1))
        pi_term = np.sum(temp1 + temp2)

        return likelihood + pi_term

    def update_alpha(self):

        def grad(alpha):
            temp = np.zeros(self.L_h)
            temp = 1/self.var_h*np.convolve(self.e_u*np.squeeze(self.h), self.g, mode='valid')/self.beta
            temp += -1/self.var_h*np.convolve(self.e_2u*np.convolve(self.g, alpha/self.beta, mode='full'), self.g, mode='valid')/self.beta

            temp += -1/(2*self.var_h)*np.convolve(self.e_2u, np.square(self.g), mode='valid')/np.square(self.beta)
            temp += (self.lamb*np.square(self.v) - alpha)*polygamma(1,alpha) - self.v/self.beta + 1

            return temp

        def hess(alpha):
            temp = np.zeros(self.L_h)
            temp = -1/self.var_h*np.convolve(self.e_2u, np.square(self.g), mode='valid')/self.beta
            temp += -polygamma(1,alpha) + (self.lamb*np.square(self.v) - alpha)*polygamma(2,alpha)

            return temp

        self.alpha = newton(grad, hess, self.alpha, 1)

def newton(grad, hess, param, num_iter):

    new_param = 0
    for i in range(num_iter):
        new_param = param - grad(param)/hess(param)
        # print('grad :', grad(param))
        # print('hess :', hess(param))
        # print('convergence :', np.mean(param - new_param))
        param = new_param

    return np.abs(param)
