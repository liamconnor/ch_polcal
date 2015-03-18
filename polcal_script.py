#############################################################################
# BETTER IDEA WHY NOT PUT ALL THE CODE INTO ONE PROGRAM AND WRITE A SCRIPT?!#
#############################################################################
import numpy as np
import h5py

import scipy.linalg.lapack

import polsol_tools as polsol
import PolcalMatrices

eig = scipy.linalg.lapack.dsyev

lr = True

PM = PolcalMatrices.PolcalMatrices([1, 0.5, 0.3, 0.0])

if lr:
     Qp = np.array([[1, -1],[1, 1]])/np.sqrt(2)
else:
     Qp = np.linalg.eig(PM.dq)[1]


f = h5py.File('outtest.hdf5', 'r')
arr = f['data'][:]

# J = np.matrix([[1., 0.01], [0.0j, 1.3]])
# J_true = J.copy()

# P_I = PM.I * PM.dic
# V_I = np.dot(J, np.dot(P_I, J.H))
# V_I = np.array(V_I)[:,:, np.newaxis]

# P_Q = np.matrix([[0.0,  PM.Q],[PM.Q, 0.0]])
# V_Q = np.array(np.dot(J, np.dot(P_Q, J.H)))[..., np.newaxis]

# P_U = np.matrix([[0.0,  1j*PM.U], [-1j*PM.U, 0.0]])
# V_U = np.array(np.dot(J, np.dot(P_U, J.H)))[..., np.newaxis]

# print P_I
# print P_Q
# print P_U

g = h5py.File('test.hdf5', 'r')

J_true = g['J'][:]
P_U = g['P_U'][:]
P_Q = g['P_Q'][:]

wp, Qp, INFO = eig(P_Q)
assert INFO==0
print "Qp"
print Qp

nfeed = 2
nfreq = arr.shape[-1]
nfreq = 1

V_I, V_Q, V_U = polsol.pop_stokes_mat(arr, nfeed=nfeed)

print 'getting unpol e-vec'
gain = polsol.get_eigenvectors(V_I, nfreq)
gain_solution = gain.copy()

#gain = polsol.normalize_evec(gain, ax=0, sqrt=False)

gain[:, :, 0] = gain[:, :, 0] / (abs(gain[:, :, 0])**2).sum(axis=0)
gain[:, :, 1] = gain[:, :, 1] / (abs(gain[:, :, 1])**2).sum(axis=0)

for freq in range(nfreq-1, nfreq):
    theta, QM = polsol.polsol_freq(gain[:, freq], V_Q[..., freq], V_U[..., freq], lr=lr)

    E = np.zeros([2, 2], np.complex128)
    E[0,0] = np.exp(1j * theta)
    E[1,1] = np.exp(-1j * theta)
    J_ = np.dot(np.dot(gain_solution[:,freq], QM), np.dot(E, np.matrix(Qp).H))
    J_ = np.matrix(J_)

    gain_solution[:, freq] = np.dot(gain_solution[:, freq], QM)

    gain_solution[:, freq, 0] = gain_solution[:, freq, 0] * np.exp(1.0J*theta)
    gain_solution[:, freq, 1] = gain_solution[:, freq, 1] * np.exp(-1.0J*theta)
    gain_solution[:, freq, :] = np.dot(gain_solution[:, freq, :], Qp)


    print "================================================================"
    print "==                 COMPARE JONES SOLUTIONS                    =="
    print "================================================================"
    print ""
    print ""

    print "True Jones matrix"
    print abs(J_true / (J_true[0,0]))
    print ""

    print "Jones solution"
    print np.round(abs(J_ / (J_[0,0])),4)
    print 

    print "================================================================"
    print "==               COMPARE RECOVERED VISIBILITIES               =="
    print "================================================================"
    print ""
    print ""


    print "True P_Q"
    print "========"
    print P_Q
    print 

    print "Recovered P_Q"
    print "============="
    P_Q_ = np.linalg.inv(J_) * \
     np.matrix(V_Q[..., freq]) * np.linalg.inv(J_.H)
    print np.round(P_Q_, 3)

    print ""
    print ""
    
    print "True P_U"
    print "========"
    print P_U

    print 

    print "RecoveredP_U"
    print "============"
    P_U_ = np.linalg.inv(J_) * \
     np.matrix(V_U[..., freq]) * np.linalg.inv(J_.H)
    print (np.round(P_U_, 4))



