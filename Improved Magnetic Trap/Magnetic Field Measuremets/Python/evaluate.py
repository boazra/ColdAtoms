# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:36:06 2017

@author: Boaz
"""



def evaluate(self, position):
    _p = np.atleast_2d(position)
    B = np.zeros(_p.shape)
    for loop in self.loops:
      B += self._evalLoop(_p, loop)
    return np.squeeze(B / self._field_units)
  
  def _evalLoop(self, p, loop):
    r_vect = (p - loop.p) * self._length_units
    r = np.linalg.norm(r_vect, axis=1, keepdims=True)
    z = r_vect.dot(loop.n.T)
    rho_vect = r_vect - np.outer(z, loop.n)
    rho = np.linalg.norm(rho_vect, axis=1)
    rho_vect[rho > self._epsilon,] = \
      (rho_vect[rho > self._epsilon,].T/rho[rho > self._epsilon]).T
    
    a = loop.r * self._length_units
    alpha2 = a*a + rho*rho + z*z - 2.*a*rho
    beta2 = a*a + rho*rho + z*z + 2.*a*rho
    beta = np.sqrt(beta2)        
    c = 4.e-7 * loop.i  # \mu_0  I / \pi
    a2b2 = alpha2 / beta2
    Ek2 = scipy.special.ellipe(1. - a2b2)
    Kk2 = scipy.special.ellipkm1(a2b2)
    
    denom = (2. * alpha2 * beta * rho)
    with np.errstate(invalid='ignore'):
      numer = c*z*((a*a + rho*rho + z*z)*Ek2 - alpha2*Kk2)
    sw = np.abs(denom) > self._epsilon
    Brho = np.zeros(numer.shape)
    Brho[sw] = numer[sw] / denom[sw]

    denom = (2. * alpha2 * beta)
    with np.errstate(invalid='ignore'):
      numer = c*((a*a - rho*rho - z*z)*Ek2 + alpha2*Kk2)
    sw = np.abs(denom) > self._epsilon
    Bz = np.full(numer.shape, np.inf)
    Bz[sw] = numer[sw] / denom[sw]

    return (Brho * rho_vect.T).T + np.outer(Bz, loop.n)
