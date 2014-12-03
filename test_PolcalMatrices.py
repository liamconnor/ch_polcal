import numpy as np
import numpy.linalg as la

import polcalqui as pq
import PolcalMatrices as P

class TestPolcalMat(P.PolcalMatrices):
	

	def test_wq_special(self):
		""" Confirm that wq has determinent 1 
			for random angles
		"""
		angles = np.random.rand(3)
		wq = self.wq(angles[0], angles[1], angles[2])

		assert np.round(la.det(wq), 5) == 1 

	def test_wq_unitary(self):
		""" Confirm that wq is unitary
		"""
		angles = np.random.rand(3)
		wq = self.wq(angles[0], angles[1], angles[2])
		dq = self.dq

		assert (np.round(wq * dq * wq.H) == dq).all()
	
	def test_mu_det(self):
		""" Check M_U's determinent
		"""
		angles = np.random.rand(3)
		mu = self.mu(angles[0], angles[1], angles[2])

		assert np.round(la.det(mu), 3) == -1

	def test_mic_det(self):
		""" Check M_IC's determinent
		"""
		angles = np.random.rand(3)
		mic = self.mic(angles[0], angles[1], angles[2])

		assert 	np.round(la.det(mic), 3) == self.I**2 - self.C**2

