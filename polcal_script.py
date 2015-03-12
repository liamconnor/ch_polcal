#############################################################################
# BETTER IDEA WHY NOT PUT ALL THE CODE INTO ONE PROGRAM AND WRITE A SCRIPT?!#
#############################################################################
import numpy as np
import h5py

import polsol_tools as polsol
import PolcalMatrices

lr = False

PM = PolcalMatrices.PolcalMatrices([1, 1, 1, 0.0])

if lr:
     Qp = np.array([[1, 1],[1, -1]])/np.sqrt(2)
else:
     Qp = np.linalg.eig(PM.dq)[1]


f = h5py.File('outtest.hdf5', 'r')
arr = f['data'][:]

nfeed = 2
nfreq = arr.shape[-1]

V_I, V_Q, V_U = polsol.pop_stokes_mat(arr, nfeed=nfeed)

"""################# SIMULATION ###############################
J = np.matrix([[ 1.+0.j,  0.53+0.05j],[ 0.8+0.j,  1.+0.j]])

print J

#J = np.matrix([[ 1.00+0.j  ,  0.10+0.3j ],[0.09+0.22j,  1.00+0.j]])
V_I = PM.I * PM.dic
V_I = np.dot(J, np.dot(V_I, J.H))
V_I = np.array(V_I)[:,:, np.newaxis]

V_Q = np.matrix([[0.0,  PM.Q],[PM.Q, 0.0]])
V_Q = np.array(np.dot(J, np.dot(V_Q, J.H)))[..., np.newaxis]

V_U = np.matrix([[0.0,  1j*PM.U],[-1j*PM.U, 0.0]])
V_U = np.array(np.dot(J, np.dot(V_U, J.H)))[..., np.newaxis]
################# SIMULATION ###############################
nfreq = V_U.shape[-1]"""

print 'getting unpol e-vec'
gain = polsol.get_eigenvectors(V_I, nfreq)
gain_solution = gain.copy()

#gain = polsol.normalize_evec(gain, ax=0, sqrt=False)

gain[:, :, 0] = gain[:, :, 0] / (abs(gain[:, :, 0])**2).sum(axis=0)
gain[:, :, 1] = gain[:, :, 1] / (abs(gain[:, :, 1])**2).sum(axis=0)

for freq in range(nfreq-1, nfreq):
    theta, QM = polsol.polsol_freq(gain[:, freq], V_Q[..., freq], V_U[..., freq], lr=lr)

    gain_solution[:, freq] = np.dot(gain_solution[:, freq], QM)

    gain_solution[:, freq, 0] = gain_solution[:, freq, 0] * np.exp(1.0J*theta)
    gain_solution[:, freq, 1] = gain_solution[:, freq, 1] * np.exp(-1.0J*theta)
    gain_solution[:, freq, :] = np.dot(gain_solution[:, freq, :], Qp)

    print (gain_solution[:, freq, :] / gain_solution[0, freq, 0])






