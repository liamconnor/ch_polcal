"""
This programs is meant to solve for the full-pol gain solution 
given two polarized calibrators (V_Q and V_U) and one unpolarized 
calibrator (V_I), which were all fitted for in Project_V.py. 

Written by Liam Connor connor@astro.utoronto.ca
"""
import numpy as np
import h5py
import cmath
import scipy.linalg.lapack

eig = scipy.linalg.lapack.dsyev

def normalize_evec(evec, ax=None, sqrt=True):
    """ Vectors are routinely normalized by their 
    length in the UL Pen polcal pipeline. Just trying 
    to emulate that.
    """

    if sqrt:
        evec[..., 0] /= np.sqrt(abs(evec[..., 0]**2).sum(axis=ax))
        evec[..., 1] /= np.sqrt(abs(evec[..., 1]**2).sum(axis=ax))
    else:
        evec[..., 0] /= abs(evec[..., 0]**2).sum(axis=ax)
        evec[..., 1] /= abs(evec[..., 1]**2).sum(axis=ax)

    return evec

def get_eigenvectors(arr, N):
    """ Step through each N freqs or times
    and diagonalize the (nfeed,nfeed) matrix.

    Parameters
    ----------
    arr : array_like
        (nfeed, nfeed, N) array
    N : np.int
        number of samples / freqs to iterate over

    """
    Gain_mat = arr.copy()

    nfeed = Gain_mat.shape[0]
                                          
    gain_approx = np.zeros([nfeed, N, 2], np.complex128)

    for k in range(N):
        # Get the last two e-vals and e-vecs. Gain_mat should be rank 2
        #w, v = np.linalg.eigh(Gain_mat[:, :, k], UPLO='U')      
        w, v, INFO = eig(Gain_mat[:, :, k])

        assert INFO==0

        v1 = v[:,-1] #/ v[0,-1]
        v2 = v[:,-2] #/ v[0,-2]

        gain_approx[:, k, -1] = cmath.sqrt(w[-1]) * v1 / np.sqrt(abs(v1**2).sum())
        gain_approx[:, k, -2] = cmath.sqrt(w[-2]) * v2 / np.sqrt(abs(v2**2).sum())

    return gain_approx

def pop_stokes_mat(arr, nfeed=2):
    """ This should take an array and create
    three Hermitian stokes matrices. Note I now 
    round at the end, since 1e-32 instead of 0.0
    seems to screw up the eigendecomposition. Who
    woulda though!

    Parameters
    ----------
    arr : array_like
        Vis array (3, ncorr, nfreq)
    nfeed : np.int
        Number of feeds

    Returns
    -------
    tuple of rounded Stokes visibilities 
    """
    V_I = np.zeros([nfeed, nfeed, arr.shape[-1]], np.complex128)
    V_Q = V_I.copy()
    V_U = V_I.copy()

    V_I[np.triu_indices(nfeed)] = arr[0]
    V_Q[np.triu_indices(nfeed)] = arr[1]
    V_U[np.triu_indices(nfeed)] = arr[2]

    V_I[np.tril_indices(nfeed)] = np.conj(arr[0])
    V_Q[np.tril_indices(nfeed)] = np.conj(arr[1])
    V_U[np.tril_indices(nfeed)] = np.conj(arr[2])

    #return np.round(V_I, 10), np.round(V_Q, 10), np.round(V_U, 10)
    return V_I, V_Q, V_U

def calc_WPW(G, V):
    """ Calculates the matrix M=WPW^\dagger using
    GVG^\dagger 
    """
    G_H = np.transpose(np.conj(G))
    M = np.dot(G_H, np.dot(V, G))

    M = np.dot(np.linalg.inv(G), np.dot(V, np.linalg.inv(G_H)))

    return M

def calc_moar(G, V_U, M_Q):
    """ This should take a gain mat and V_Q for each
    frequency.

    """
    #wM, QM = np.linalg.eig(M_Q)
    wM, QM, INFO = eig(M_Q)

    assert INFO==0
    

    QM = normalize_evec(QM, sqrt=True)
    QSU = np.dot(G, QM)
    QSU_H = np.transpose(np.conj(QSU))

    # Really do need to double check all of this
    WP3W = np.dot(np.dot(QSU_H, V_U), QSU)

    return WP3W, QM

def solve_phase(WP3W, lr=False):
    phase = WP3W[0, 1] # Assumes the xy basis

    if lr:
        phase = 1.0J*WP3W[0,1]

    theta = np.angle(phase)/2. # Should solve equation 4.13 in Greg's thesis
    
    return theta

def polsol_freq(gain, V_Q, V_U, lr=False):
    """ Needs to return QM
    """
    M_Q = calc_WPW(gain, V_Q)
    WP3W, QM = calc_moar(gain, V_U, M_Q)
    theta = solve_phase(WP3W, lr=lr)
    
    
    print "theta =", theta

    return theta, QM


"""
====================================
= Couple notes on the pol matrices =
====================================

1) Ginv V_i Ginv.H = W.H P_i W \equiv M_i
2) M_(QU)[0,0] = -M_(QU)[1,1]
3) M_Q[0,1] = M_Q[1,0]
4) M_Q[0,0] = Q*cos(2 alph)
5) M_Q[0,1] = Q*cos(alph)*sin(alph)*e^{i(beta-gamma)}
6) M_U[0,1] = conj(M[1,0])



"""




