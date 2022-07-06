#------------------------------------------------------------------------------
# Calculate tidal Love numbers for model planets
#------------------------------------------------------------------------------
import numpy as np

def lovek2(zvec, dvec):
    """ Tidal Love number k2 from density profile.

    Inputs
    -------
    zvec : real finite positive vector
      Normalized radii starting near (not at) center
    dvec : real finite nonnegative vector
      Density at corresponding radius

    Outputs
    --------
    k2 : scalar, real
      Love number k2

    Algorithm
    ---------
    Sterne 1939 as described in Buhler 2016.
    """
    assert zvec.ndim == dvec.ndim == 1
    assert zvec.size == dvec.size
    assert all(zvec > 0)
    assert all(dvec >= 0)

    # Climb up the radius with Buhler (2016) eq. 2
    m = dvec[0]*zvec[0]**3; # starting mass
    rhom = m/zvec[0]**3; # starting mean density
    eta = 0
    for k in range(zvec.size - 1):
        s1 = (6 - 6*(dvec[k]/rhom)*(eta + 1) + eta - eta**2)/zvec[k];
        zhalf = zvec[k] + 0.5*(zvec[k+1] - zvec[k]);
        dhalf = dvec[k] + 0.5*(dvec[k+1] - dvec[k]);
        mhalf = m + dhalf*(zhalf**3 - zvec[k]**3);
        rhalf = mhalf/zhalf**3;
        ehalf = eta + s1*(zhalf - zvec[k]);
        s2 = (6 - 6*(dhalf/rhalf)*(ehalf + 1) + ehalf - ehalf**2)/zhalf;
        eta = eta + s2*(zvec[k+1] - zvec[k]);
        m = mhalf + dvec[k+1]*(zvec[k+1]**3 - zhalf**3);
        rhom = m/zvec[k+1]**3;
    k2 = (3 - eta)/(2 + eta)
    return k2

if __name__ == '__main__':
    # validate on uniform density
    N = 2048
    eps = np.spacing(1)
    zvec = np.linspace(1/N, 1, N)
    dvec = np.ones_like(zvec)
    k2_calc = lovek2(zvec,dvec)
    k2_expect = 1.5
    k2_err = abs(k2_calc - k2_expect)/abs(k2_expect + eps)
    print(f"Uniform density, N={N}")
    print(f"expected: {k2_expect:g}; calculated: {k2_calc:g}; err: {k2_err:g}")

    # validate with n=1 polytrope approximating Jupiter-like planet
    pi = np.pi
    M = 317.8*5.9722e+24
    R = 71492e3
    G = 6.67430e-11
    K = 2*G/pi*R**2
    a = np.sqrt(2*pi*G/K)
    R = pi/a; # should reproduce R
    rho_av = 3*M/(4*pi*R**3)
    rho_c = (pi**2/3)*rho_av
    r = np.linspace(1/N,1,N)*R
    rho_exact = rho_c*np.sin(a*r)/(a*r)
    rho_exact[-1] = 0 # could end up on the negative side of double zero...
    k2_calc = lovek2(r,rho_exact)
    k2_expect = (15/pi**2) - 1
    k2_err = abs(k2_calc - k2_expect)/abs(k2_expect + eps);
    print()
    print(f"n=1 polytrope, N={N}")
    print(f"expected: {k2_expect:g}; calculated: {k2_calc:g}; err: {k2_err:g}")
