"""
This programs is meant to solve for the full-pol gain solution given two polarized calibrators (V_Q and V_U) and one unpolarized calibrator (V_I), which were all fitted for in project_V.py. 
written by Liam Connor connor@astro.utoronto.ca
"""
basis='lr'

import numpy as np
import h5py
import cmath

nfreq = 1000
n_ant = 2
ncorr = n_ant*(n_ant+1)/2

#f = np.fromfile('dataPOL.dat', np.complex64)
#f = f.reshape(-1, ncorr, nfreq)

fil = h5py.File('test.hdf5','r')
f = fil['arr'][:]

print f.shape

V_I = f[0]
V_Q = f[1]
V_U = f[2]

if basis is 'lr':
	Qp = np.array([[-1,1],[1,1]])/np.sqrt(2)
if basis is 'xy':
	Qp = np.array([[1,0],[0,1]])

def antmap(i,j,n):
        i=i-1;j=j-1
        if i>j:
                return ((2*n*j - j**2 + j) / 2) + (i-j)
        else:
                return ((2*n*i - i**2 + i) / 2) + (j-i)

def get_eigenvectors(Data, N):
        """                                                                                                                        Function to get eigenvectors in LL and RR for each time or frequency                                                                                                                                                                                  ====params====                                                                                                                                 
        Data - Array with unpol visibilites 
        N - Either number of frequencies or number of times on which to solve e-vec                                                                    
        """
        Gain_mat = np.zeros([n_ant, n_ant, N], np.complex64)

        for ant_i in range(n_ant):
                for ant_j in range(n_ant):
                        if ant_i <= ant_j:
                                Gain_mat[ant_i,ant_j,:] = Data[antmap(ant_i+1, ant_j+1, n_ant), :]
                                Gain_mat[ant_j,ant_i,:] = np.conj(Data[antmap(ant_i+1, ant_j+1, n_ant), :]) #Should be Hermitian                  
                                Gain_mat[ant_j,ant_j,:] = 0 # Set autocorrelations to zero            

        # Now try to find eigenvectors of 64x64 gainmatrices for all N
                                              
        gain_approx = np.zeros([n_ant, N, 2], np.complex64)

        for k in range(N):
                # Get the last two e-vals and e-vecs. Gain_mat should be rank 2
                w, v = np.linalg.eigh(Gain_mat[:, :, k], UPLO='U')		

                v1 = v[:,-1] / v[0,-1]
                v2 = v[:,-2] / v[0,-2]

                gain_approx[:, k, 0] = cmath.sqrt(w[-1]) * v1 / np.sqrt(abs(v1**2).sum())
                gain_approx[:, k, 1] = cmath.sqrt(w[-2]) * v2 / np.sqrt(abs(v2**2).sum())

        return gain_approx

print 'getting unpol e-vec'
gain = get_eigenvectors(V_I, nfreq)
gain_solution = get_eigenvectors(V_I, nfreq)

gain[:,:,0] = gain[:,:,0] / ((abs(gain[:,:,0])**2).sum(axis=0))
gain[:,:,1] = gain[:,:,1] / ((abs(gain[:,:,1])**2).sum(axis=0))

print 'beginning pol part'
for freq in range(nfreq):
#	print freq
	QM = np.zeros([2, 2],np.complex64)

	V_Q_ = V_Q[:,freq]
	V_U_ = V_U[:,freq]
	V_Qmat = np.zeros([n_ant, n_ant], np.complex64)
	V_Umat = np.zeros([n_ant, n_ant], np.complex64)

        for ant_i in range(n_ant):
		for ant_j in range(n_ant):
			if ant_i < ant_j:
				V_Qmat[ant_i, ant_j] = V_Q_[antmap(ant_i+1, ant_j+1, n_ant)]
				V_Qmat[ant_j, ant_i] = np.conj(V_Q_[antmap(ant_i+1, ant_j+1, n_ant)])

				V_Umat[ant_i, ant_j] = V_U_[antmap(ant_i+1, ant_j+1, n_ant)]
                                V_Umat[ant_j, ant_i] = np.conj(V_U_[antmap(ant_i+1, ant_j+1, n_ant)])

                             	V_Qmat[ant_j,ant_j] = 0 # Set autocorrelations to zero
				V_Umat[ant_j,ant_j] = 0

	M = np.dot(np.dot(np.transpose(np.conj(gain[:, freq, :])),V_Qmat),gain[:, freq, :]) # This is WPW in Ue-Li's code
	
	wM, QM = np.linalg.eigh(M,UPLO='U') #QM equiv to wt in Ue-Li's code
	
	QM[:,0] = QM[:,0] / np.sqrt(abs(QM[:,0]**2).sum())
	QM[:,1]= QM[:,1] / np.sqrt(abs(QM[:,1]**2).sum())

	QSU = np.dot(gain[:,freq,:], QM)
	gain_solution[:, freq, :] = np.dot(gain_solution[:, freq, :], QM)

	WP3W = np.dot(np.dot(np.transpose(np.conj(QSU)), V_Umat), QSU) # Should be equal to [[0,-iUexp(2i\theta)][iUexp(-2i\theta,0]]
	if basis is 'lr':
		phase = -1.0J*WP3W[0,1] 
        if basis is 'xy':
		phase = WP3W[0,1]
	
	theta = np.angle(phase)/2. # Should solve equation 4.13 in Greg's thesis
	print "theta =",theta

	gain_solution[:, freq, 0] = gain_solution[:, freq, 0] * np.exp(1.0J*theta)
	gain_solution[:, freq, 1] = gain_solution[:, freq, 1] * np.exp(-1.0J*theta)
	gain_solution[:, freq, :] = np.dot(gain_solution[:, freq, :], Qp)
	
	print gain_solution[0, freq, :]

gain_solution.tofile('gain2.liam.dat')




	
	

