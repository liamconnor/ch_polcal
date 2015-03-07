import numpy as np
import h5py

import ch_util.ephemeris as eph
import coord_tools 

# What does this program need to do?
# it needs to get the parallactic angle as a fct of time
# Needs to read in a (nfreq, ncorr, ntimes) array
# it needs to generate A and B for Ax = B
# it needs to solve for x

# it needs to save it all down
AROLONGITUDE = -78.0730
AROLATITUDE = 45.9555

RA_src, dec_src = 53.51337, 54.6248916

def drao_aro_RA(times):
    """ Go from DRAO transit RA to ARO transit RA
    """
    RA_trans = eph.transit_RA(times)
    long_diff = AROLONGITUDE - eph.CHIMELONGITUDE 
    
    return RA_trans + long_diff

def ha_src(times, RA_src):
    """ Given a unix time and a source RA
    get its hour angle
    """
    RA_trans = drao_aro_RA(times)
    return RA_trans - RA_src

def get_parallactic(times, RA_src, dec_src):
    """ Takes unix times, source location
    and return parallactic angle in degrees
    """
    ha_r = np.radians(ha_src(times, RA_src))
    lat_r = np.radians(AROLATITUDE)
    dec_r = np.radians(dec_src)

    par_r = coord_tools.local_coords_dhl(dec_r, ha_r, lat_r)[-1]

    return np.degrees(par_r)


if __name__=='__main__':
    fname = 'test.hdf5'

    f = h5py.File(fname, 'r')
    data = f['arr'][:]
    times = f['times'][:]
    
    rmodel = abs(data.copy())

    nfreq = data.shape[0]
    ncorr = data.shape[1]
    ntimes = data.shape[-1]
    
    assert len(times)==ntimes
    
    phase = 2 * get_parallactic(times, RA_src, dec_src)
    phase = -2 * coord_tools.local_coords_dhl\
        (np.radians(dec_src), times, np.radians(AROLATITUDE))[-1]

    A = np.zeros([ntimes, 3], np.complex128)
    x_sol = np.zeros([3, ncorr, nfreq], np.complex128)

    for corr in range(ncorr):
        A[:, 0] = rmodel[:, corr, nfreq/2] / rmodel[:, corr, nfreq/2].sum()
        A[:, 1] = np.exp(1.0J * phase) * rmodel[:, corr, nfreq/2] / rmodel[:, corr, nfreq/2].sum()
        A[:, 2] = np.exp(-1.0J * phase) * rmodel[:, corr, nfreq/2] / rmodel[:, corr, nfreq/2].sum()

        B = data[:, corr, :]

        x_sol[:, corr,:] = np.linalg.lstsq(A, B)[0]
        
        print corr, np.round(x_sol[1, corr] + x_sol[2, corr], 3)
        print 1j * np.round(x_sol[1, corr] - x_sol[2, corr], 3)
        print x_sol[:, corr, 0]
        print x_sol.shape

    dataPOL = np.zeros_like(x_sol)

    dataPOL[0,:,:] = x_sol[0,:,:] 
    # Should be unpolarized solution for each freq and baseline                   
    dataPOL[1,:,:] = x_sol[1,:,:] + x_sol[2,:,:]
    # Gives us Stokes Q visibities                                                                             
    dataPOL[2,:,:] = 1.0J*(x_sol[1,:,:] - x_sol[2,:,:]) 
    # Gives us Stokes U    
    
#    print (dataPOL[0]), dataPOL[1], (dataPOL[2])

