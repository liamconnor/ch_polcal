import numpy as np


class PolcalMatrices:

	def __init__(self, stokes):
		self.I, self.Q, self.U, self.C = stokes

	@property 
	def qp(self):
		return np.matrix([[1+1j, 1j-1],[1+1j, 1-1j]]) / 2.0

	@property 	
	def dq(self):
		return np.matrix([[-1, 0], [0, 1]]) # Should this have det -Q**2 or -1?

	@property
	def du(self):
		return np.matrix([[0, 1], [1, 0]])  # Should this have det -U**2 or -1? 

	@property	
	def dic(self):
		return np.matrix([[self.I, 1j*self.C], [-1j*self.C, self.I]])
	
	def polmat(self, di):
		return np.dot(self.qp, np.dot(di, self.qp.H))

	def phom2ab(self, alpha, beta):
		""" Reparameterizes phase angles

		Parameters
		----------
		alpha : float, array_like
			angle, alpha
		beta : float, array_like
			angle, beta

		Returns
		-------
		ph : float, array_like
			angle, phi
		om : float, array_like
			angle, omega
		"""
		ph = (alpha + beta) / 2.
		om = (alpha - beta) / 2.

		return ph, om

	def wq(self, ph, om, th, matrix=True):
		""" Generate an SU(1,1) matrix 

		Parameters  
		----------
		ph : float, array_like
			angle, phi
		om : float, array_like
			angle, omega
		th : float, array_like
			angle, theta

		Returns 
		-------
		wq : array_like, matrix
			SU(1,1) matrix
		"""

		if matrix is True:
			a = np.exp(1j*ph) * np.cosh(th)
			b = np.exp(1j*om) * np.sinh(th)
			
			return np.matrix([[a, b], [np.conj(b), np.conj(a)]])

		elif matrix is False:
			ph = np.array(ph)
			if ph.shape != ():
				ph = ph[:, np.newaxis, np.newaxis]
			th = np.array(th)
			if th.shape != ():
				th = th[np.newaxis, :, np.newaxis]
			
			om = np.array(om)
			if om.shape != ():
				th = om[np.newaxis, np.newaxis, :]

			a = np.exp(1j*ph) * np.cosh(th)
			b = np.exp(1j*om) * np.sinh(th)
			
			return np.array([[a, b], [np.conj(b), np.conj(a)]])
		else:
			raise Exception("What kind of array do you want returned?")

	def mu(self, ph, om, th, matrix=True):
		""" Generate M_U matrix 

		Parameters  
		----------
		ph : float, array_like
			angle, phi
		om : float, array_like
			angle, omega
		th : float, array_like
			angle, theta
		matrix : bool, optional
			method will return a matrix if this is True

		Returns 
		-------
		mu : array_like, matrix
			M_U = G**(-1)*V_U*G**(-\dagger) 
		"""
		
		a = np.cos(ph - om) * np.sinh(2*th)
		b = np.exp(2j*ph)*np.cosh(th)**2 + np.exp(2j*om)*np.sinh(th)**2

		if matrix is True:
			return np.matrix([[a, b], [np.conj(b), a]])

		elif matrix is False:
			return np.array([[a, b], [np.conj(b), a]])

		else:
			raise Exception("What kind of array do you want returned?")

	def mic(self, ph, om, th):
		""" Generate M_IC matrix 

		Parameters  
		----------
		ph : float, array_like
			angle, phi
		om : float, array_like
			angle, omega
		th : float, array_like
			angle, theta
		matrix : bool, optional
			method will return a matrix if this is True

		Returns 
		-------
		mic : array_like, matrix
			M_IC = G**(-1) * V_{IC} * G**(-\dagger) 
		"""
		a = np.cosh(th) * np.exp(1j*ph)
		b = np.sinh(th) * np.exp(1j*om)

		m00 = (abs(a)**2 + abs(b)**2)*self.I + 1j*self.C*(a*np.conj(b) - np.conj(a)*b)
		m01 = a*b*self.I - 1j*b**2 * self.C + 1j*a**2*self.C + self.I*a*b

		return np.matrix([[m00, m01],[np.conj(m01), m00]])
