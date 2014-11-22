import numpy as np
import polcalqui as pq
######


#######============

def A02(M, p_o): #checked
	phio, omo, tho = p_o

	return 2 * np.sinh(2*tho) * M[0,0] - 2 * (np.exp(1j*(phio + omo)) * M[1,0] * np.cos(2*tho)).real

def A00(M, p_o): #
	phio, omo, tho = p_o

	return -2 * (  1j * np.exp(1j*(phio + omo)) * M[1,0] * np.sinh(tho) * np.cosh(tho) ).real

def A01(M, p_o):
	phio, omo, tho = p_o

	return -2 * (  1j * np.exp(1j*(phio + omo)) * M[1,0] * np.sinh(tho) * np.cosh(tho) ).real

#######==============

def A20(M, p_o):
	phio, omo, tho = p_o

	return 2j*np.cosh(tho)*np.sinh(tho)*np.exp(1j*(omo-phio))*M[0,0] \
	- 2j*np.cosh(tho)**2*np.exp(-2j*phio)*M[0,1]

def A21(M, p_o):
	phio, omo, tho = p_o

	return -2 * np.exp(1j*(omo - phio)) * 1j * np.sinh(tho)*np.cosh(tho)*M[0,0] \
	+ np.exp(2j*omo)*2j*np.sinh(tho)**2*M[1,0]

def A22(M, p_o):
	phio, omo, tho = p_o

	return -2 * np.exp(1j*(omo-phio))*np.cosh(2*tho)*M[0,0] \
	+ np.exp(2j*omo)*np.sinh(2*tho)*M[1,0] + np.exp(-2j*phio)*np.sinh(2*tho)*M[0,1]

#######==============

def y0(M, p_o): # Double checked, seems good.
	phio, omo, tho = p_o

	a = M[0,0] * (np.cosh(2*tho) - 2*tho*np.sinh(2*tho))
	b = -2 * (np.exp(1j*(phio + omo)) * M[1,0] * \
		(np.sinh(tho) * np.cosh(tho) - tho * np.cos(2*tho) - \
			1j*np.sinh(tho)*np.cosh(tho)*(phio + omo))).real
	return -(a + b)

def y1(M, p_o):
	phio, omo, tho = p_o

	a = -2*np.exp(1j*(omo - phio)) * (np.sinh(tho)*np.cosh(tho) + 1j*(phio - omo)*np.sinh(tho)*np.cosh(tho) - tho*np.cosh(2*tho))*M[0,0]
	b = np.exp(2j*omo)*( np.sinh(tho)**2 - tho*np.sinh(2*tho) - 2j*np.sinh(tho)**2*omo )*M[1,0]
	c = np.exp(-2j*phio)*( np.cosh(tho)**2+2j*phio*np.cosh(tho)**2 - tho*np.sinh(2*tho))*M[0,1]

	return -(a+b+c)

######=================


def populate_mat(p, p_o):
	A = np.zeros([3,3], np.complex128)

	MU = pq.MU1(p[0], p[1], p[2])
	MIV = pq.MIV(p[0], p[1], p[2], 1, 0.5)

	A[0,0] = (A00(MU, p_o)).real
	A[0,1] = (A01(MU, p_o)).real
	A[0,2] = (A02(MU, p_o)).real

	A[1,0] = (A20(MU, p_o)).imag	
	A[1,1] = (A21(MU, p_o)).imag
	A[1,2] = (A22(MU, p_o)).imag

	A[2,0] = (A20(MIV, p_o)).real
	A[2,1] = (A21(MIV, p_o)).real
	A[2,2] = (A22(MIV, p_o)).real
	
	RHS = y0(MU, p_o).real, y1(MU, p_o).imag, y1(MIV, p_o).real
	x = np.matrix(p_o).T 
	A = np.matrix(A)

	return A, A*x, np.matrix(RHS).T

def check(A, p):
	x = np.matrix(p).T 
	A = np.matrix(A)

	return A*x

