import numpy as np
import unittest
import time
import h5py

import Project_vis as pv
import ch_util.ephemeris as eph
import PolcalMatrices as PM
import coord_tools

class Test_Project_vis(unittest.TestCase):

    t0 = time.time()
    ntimes = 1000
    times = np.linspace(t0-3600, t0+3600, ntimes)
    RA_src, dec_src = 53.51337, 54.6248916

    def test_drao_aro_RA(self):
        times = self.times
        diff_true = np.round(np.mod(eph.CHIMELONGITUDE - (-78.0730), 360))

        DRAO_RA_trans = eph.transit_RA(times)
        ARO_RA_trans = pv.drao_aro_RA(times)

        diff = np.round(np.mod(DRAO_RA_trans - ARO_RA_trans, 360).mean())

        assert diff == diff_true
        assert np.round(np.std(diff), 2)==0.0

    def test_ha_src(self):
        times = self.times
        ha_src = pv.ha_src(times, self.RA_src)
        
        diff_true = np.round(np.mod(eph.CHIMELONGITUDE - (-78.0730), 360))
        trans_time = eph.transit_times(self.RA_src + diff_true, times[0])
        ha_src_trans = np.mod(np.round(pv.ha_src(trans_time, self.RA_src)), 360)

        ha_src_pretrans = pv.ha_src(trans_time - 1000, self.RA_src)

        assert ha_src_trans == 0.0
        assert (ha_src < 360).all()
        assert ha_src_pretrans < 0.0
    
    def test_unix_to_ha(self):

        pass

    def test_get_parallactic(self):
        times = self.times
        
        par = pv.get_parallactic(times, self.RA_src, self.dec_src)
        
        assert len(par)==len(times)
        assert (par < 360.0).all()
        
        

def _test_simulate_observation():
    """ Simulate a time stream and get polcal solution
    """
    lr = False
    I, Q, U, V = 1, 0.05, 0.15, 0.0
    PolMat = PM.PolcalMatrices([I, Q, U, V])

    RA_src, dec_src = 53.51337, 54.6248916

    ntimes = 1000

    ha_r = np.linspace(-np.pi/4, np.pi/4, ntimes)
    phase = 2 * coord_tools.local_coords_dhl(np.radians(dec_src), \
               ha_r, np.radians(pv.AROLATITUDE))[-1][:, np.newaxis, np.newaxis]

    J = np.array([[1., 0.3+.1j], [0.2+0.18j, 1.2]]) #* np.exp(1j*33*1.0)
    print J / J[0,0]

    P_Q = np.array(Q * PolMat.dq).astype(np.complex128)
    P_U = np.array(U * PolMat.du).astype(np.complex128)
    P_I = np.array(I * PolMat.dic).astype(np.complex128)

    if lr:
        P_U[1,0] *= -1.0j
        P_U[0,1] *= 1.0j
        P_Q *= 0.0
        P_Q[0,1] = Q 
        P_Q[1,0] = Q

        print "P_U, P_Q" 
        print P_U
        print P_Q 
        print ""

    V_Q = np.dot(J, np.dot(P_Q, np.conj(J.transpose())))[np.newaxis]
    V_U = np.dot(J, np.dot(P_U, np.conj(J.transpose())))[np.newaxis]
    V_I = np.dot(J, np.dot(P_I, np.conj(J.transpose())))[np.newaxis]

    print ""
    print "I : xx, xy, yy", np.round([V_I[0,0,0], V_I[0,0,1], V_I[0,1,1]], 3)
    print "Q : xx, xy, yy", np.round([V_Q[0,0,0], V_Q[0,0,1], V_Q[0,1,1]], 3)
    print "U : xx, xy, yy", np.round([V_U[0,0,0], V_U[0,0,1], V_U[0,1,1]], 3)
    print ""

    noise = np.random.normal(0, 0.1, ntimes*2*2).reshape(-1, 2, 2)\
        + 1j * np.random.normal(0, 0.1, ntimes*2*2).reshape(-1, 2, 2)

    V = V_I + 0.5 * (V_Q - 1.0J*V_U) *  np.exp(1.0J*phase) \
            + 0.5*(V_Q + 1.0J*V_U)*np.exp(-1.0J*phase) + 0. * noise

    assert np.round(np.linalg.det(P_Q), 2) == np.round(-1 * Q**2, 2)
    assert np.linalg.det(P_U).astype(np.float) == -1 * U**2
    assert np.linalg.det(P_I).astype(np.float) == 1 * I**2
        
    return V.reshape(-1, 4)[:, (0,1,3)], ha_r, phase
        
if __name__=='__main__':
    V, ha_r, ph = _test_simulate_observation()
    V = np.array(V[..., np.newaxis] * np.ones([1, 1, len(ha_r)]))
    
    f = h5py.File('test.hdf5','w')
    f.create_dataset('arr', data=V)
    f.create_dataset('times', data=ha_r)
    f.create_dataset('phase', data=ph[:, 0, 0])
    
    f.close()

    unittest.main()
    
    
