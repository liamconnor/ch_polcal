import numpy as np

def Wq1(ph, om, th, matrix='on'):
    a = np.exp(1j*ph) * np.cosh(th)
    b = np.exp(1j*om) * np.sinh(th)

    if matrix=='on':
        return np.matrix([[a, b], [np.conj(b), np.conj(a)]])

    elif matrix=='off':
        ph = np.array(ph)
        if ph.shape!=():ph = ph[:, np.newaxis, np.newaxis]
        th = np.array(th)
        if th.shape!=():th = th[np.newaxis, :, np.newaxis]
        om = np.array(om)
        if om.shape!=():th = om[np.newaxis, np.newaxis, :]
        a = np.exp(1j*ph) * np.cosh(th)
        b = np.exp(1j*om) * np.sinh(th)
        print a.shape, b.shape
        return np.array([[a, b], [np.conj(b), np.conj(a)]])
    else:
        raise Exception("WHAT KIND OF ARRAY DO YOU WANT?")


def Wq(alpha, beta, th, matrix='on'):
    a = np.exp(1j * (alpha + beta) / 2.) * np.cosh(th)
    b = np.exp(1j * (alpha - beta) / 2.) * np.sinh(th)

    if matrix=='on':
        return np.matrix([[a, b], [np.conj(b), np.conj(a)]])

    elif matrix=='off':
        alpha = np.array(alpha)
        if alpha.shape!=():alpha = alpha[:, np.newaxis, np.newaxis]
        th = np.array(th)
        if th.shape!=():th = th[np.newaxis, :, np.newaxis]
        beta = np.array(beta)
        if beta.shape!=():beta = beta[np.newaxis, np.newaxis, :]
        a = np.exp(1j * (alpha + beta) / 2.) * np.cosh(th)
        b = np.exp(1j * (alpha - beta) / 2.) * np.sinh(th)
        print a.shape, b.shape
        return np.array([[a, b], [np.conj(b), np.conj(a)]])
    else:
        raise Exception("WHAT KIND OF ARRAY DO YOU WANT?")

def DD00(M, p_o, p):
    """
    Parameters
    ==========
    M: complex matrix
       Should be the matrix M_i = Ginv * D_i * GHinv
    p_o: 
       angles we expand around. Ordered [phi, omega, theta]
    p:
       [phi, omega, theta]
    
    Returns
    =======
    The 01 element of the product matrix: la.inv(W) * M * la.inv(W.H)
    """
    phio, omo, tho = p_o
    phi, om, th = p

    a = (np.cosh(2*tho) + 2*(th - tho) * np.sinh(2*tho)) * M[0,0]

    b = np.exp(1j * (phio + omo)) 

    c = (np.sinh(tho) * np.cosh(tho) + (th - tho) * np.cos(2*tho)) + np.sinh(tho) * np.cosh(tho) * (1j * (phi + om - phio - omo))
    
    return a - 2 * (b * c * M[1,0]).real

def ao(M, p_o, p):
    phio, omo, tho = p_o
    phi, om, th = p

    return M[0,0]*2*np.sinh(2*tho)

def a1(M, p_o, p):
    phio, omo, tho = p_o
    phi, om, th = p
    
    c = (np.sinh(np.sinh(th) * np.cosh(tho) * (1j)))
    
    return (-2 * np.real(np.exp(1j*(phio + omo))*( np.sinh  )))


def DD01(M, p_o, p):
    """
    Parameters
    ==========
    M: complex matrix
       Should be the matrix M_i = Ginv * D_i * GHinv
    p_o: 
       angles we expand around. Ordered [phi, omega, theta]
    p:
       [phi, omega, theta]
    
    Returns
    =======
    The first element of the product matrix: la.inv(W) * M * la.inv(W.H)
    """
    
    phio, omo, tho = p_o
    phi, om, th = p
    a = -2*np.exp(1j*(omo - phio)) * ( np.sinh(tho)*np.cosh(tho) + 1j*(om-phi-omo+phio)*np.sinh(tho)*np.cosh(tho) + (th-tho)*np.cosh(2*tho) ) * M[0,0]
    b = np.exp(2j*omo)*(np.sinh(tho)**2 + (th - tho)*np.sinh(2*tho) + 2j*np.sinh(tho)**2*(om-omo))*M[1,0]
    c = np.exp(-2j*phio)*(np.cosh(tho)**2 - 2j*(phi-phio)*np.cosh(tho)**2 + (th-tho)*np.sin(2*tho))*M[0,1]
    
    return a + b + c

#    a = -2 * ( np.exp(1j*(omo - phio)*(1 + 1j*(om - phi - omo + phio)))) * M[0,0]
#    b = np.exp(-2j*phio) * (np.cosh(tho)**2 + 2*np.cosh(tho)*np.sinh(tho)*(th - tho) - 2j*np.cosh(tho)**2*(phi-phio)) * M[0,1]
#    c = np.exp(2j*omo) * (np.sinh(tho)**2 + 2*np.cosh(tho)*np.sinh(tho)*(th - tho) + 2j*np.sinh(tho)**2*(om - omo)) * M[1,0]


def MU1(ph, om, th, matrix='on'):
    a = np.cos(ph - om) * np.sinh(2*th)
    b = np.exp(2j*ph)*np.cosh(th)**2 + np.exp(2j*om)*np.sinh(th)**2
    if matrix=='on':
        return np.matrix([[a, b],[np.conj(b), a]])
    elif matrix=='off':
        return np.array([[a, b],[np.conj(b), a]])
    else:
        raise Exception("WHAT KIND OF ARRAY DO YOU WANT?")

def MU(alpha, beta, theta):
    alpha = alpha[:, np.newaxis, np.newaxis]
    beta = beta[np.newaxis, :, np.newaxis]
    theta = theta[np.newaxis, np.newaxis]

    a = np.cos(beta) * np.sinh(2*theta)
    b = np.cosh(theta)**2 * np.exp(1j*(alpha + beta)) + np.sinh(theta)**2 * np.exp(1j*(alpha - beta))
    
    return np.array([[a, b],[np.conj(b), a]])

def MI(th, alpha, be, A, V, matrix='on'):
    a = A*np.cosh(2*th) - 2*V*np.sin(be)*np.sinh(2*th)
    b = 1j*V*(np.exp(1j*(alpha + be))*np.cosh(th)**2 - np.exp(1j*(alpha - be))*np.sinh(th)**2 + A*np.exp(1j*(alpha))*np.sinh(2*th))

    if matrix=='on':
        return np.matrix([[a, b],[np.conj(b), a]])
    elif matrix=='off':
        ph = np.array(ph)
        if ph.shape!=():ph = ph[:, np.newaxis, np.newaxis]
        th = np.array(th)
        if th.shape!=():th = th[np.newaxis, :, np.newaxis]
        om = np.array(om)
        if om.shape!=():th = om[np.newaxis, np.newaxis, :]
        a = A*np.cosh(2*th) - 2*V*np.sin(be)*np.sinh(2*th)
        b = 1j*V*(np.exp(1j*(alpha + be))*np.cosh(th)**2 - np.exp(1j*(alpha - be))*np.sinh(th)**2 + A*np.exp(1j*(alpha))*np.sinh(2*th))
        return np.array([[a, b],[np.conj(b), a]])
    else:
        raise Exception("WHAT KIND OF ARRAY DO YOU WANT?")

def MIV(ph, om, th, A, V):
    a = np.cosh(th)*np.exp(1j*ph)
    b = np.sinh(th)*np.exp(1j*om)

    M00 = (abs(a)**2 + abs(b)**2)*A + 1j*V*(a*np.conj(b) - np.conj(a)*b)
    M01 = a*b*A - 1j*b**2 * V + 1j*a**2*V + A*a*b

    return np.matrix([[M00, M01],[np.conj(M01), M00]])

def M01_(M_I, alpha, beta, theta):
    a = np.exp(1j * (alpha[:,np.newaxis,np.newaxis] + beta[np.newaxis, :, np.newaxis]) / 2.) * np.cosh(theta[np.newaxis, np.newaxis])
    b = np.exp(1j * (alpha[:,np.newaxis,np.newaxis] - beta[np.newaxis, :, np.newaxis]) / 2.) * np.sinh(theta[np.newaxis, np.newaxis])
    return 2 * a * b * M_I[0,0] + b**2 * M_I[1,0] + a**2 * M_I[0,1]

"""
def MI_antimation():
    fig1 = plt.figure()
    plt.xlabel('beta [pi]')
    plt.ylabel('alpha [pi]')
    theta = np.linspace(0, 2*np.pi, 100)
    images = []
    for th in theta:
        images.append((plt.imshow(M_[:, :, th].real, extent=[0, 2, 0, 2], vmin=0.1*(M_[:, :, th].real).min(), vmax=10*(M_[:, :, th].real).min()))
    
    im_ani = animation.ArtistAnimation(fig1, images, interval=50, repeat_delay=3000, blit=True)
    plt.show()

th, bet, alph = 2*np.pi*np.random.rand(3)
A, V = 1.0, 0.10

M_I = MI(th, bet, alph, A, V)
M_ = M01_(M_I, alpha, beta, theta)

fig1 = plt.figure()
plt.xlabel('beta [pi]')
plt.ylabel('alpha [pi]')
theta = np.linspace(0, 2*np.pi, 200)
images = []

for th in theta:
    images.append((plt.imshow(abs(np.angle(M_[:, :, th])), extent=[0, 2, 0, 2], vmin=np.pi/2. - .2, vmax=np.pi/2 + .2),))
#    images.append((plt.imshow(abs((M_[:, :, th].real - (M_[:, :, th].real).min())), extent=[0, 2, 0, 2], vmin=0, vmax=1000),))

plt.colorbar()
im_ani = animation.ArtistAnimation(fig1, images, interval=50, repeat_delay=3000, blit=True)

plt.show()
"""
