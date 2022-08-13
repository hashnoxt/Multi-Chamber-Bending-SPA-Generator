
import numpy as np


# --------------------------------------------------------------------------------
# Material model name.
# --------------------------------------------------------------------------------
def name(): return 'ogden3'
def pname(): return 'Ogden-3'
def params(): return 'u1 a1 u2 a2 u3 a3'
def descr(): return 'Ogden Model with order 3 (modified form).'


# NOTE - this is the Abaqus form of the functions. Ogden2004 is similar, but they
#        show these functions as being multiplied by (a[i]/2.0)

# --------------------------------------------------------------------------------
# Function defining the uniaxial stress given strain.
# --------------------------------------------------------------------------------
def stressU(x, u1, a1, u2, a2, u3, a3):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L, a1-1.0) - np.power(L, -0.5*a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L, a2-1.0) - np.power(L, -0.5*a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L, a3-1.0) - np.power(L, -0.5*a3-1.0)) / a3
    return S1+S2+S3


# --------------------------------------------------------------------------------
# Function defining the biaxial stress given strain.
# --------------------------------------------------------------------------------
def stressB(x, u1, a1, u2, a2, u3, a3):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L, a1-1.0) - np.power(L, -2.0*a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L, a2-1.0) - np.power(L, -2.0*a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L, a3-1.0) - np.power(L, -2.0*a3-1.0)) / a3
    return S1+S2+S3


# --------------------------------------------------------------------------------
# Function defining the planar stress given strain.
# --------------------------------------------------------------------------------
def stressP(x, u1, a1, u2, a2, u3, a3):
    L = 1.0+x
    S1 = 2.0 * u1 * (np.power(L, a1-1.0) - np.power(L, -a1-1.0)) / a1
    S2 = 2.0 * u2 * (np.power(L, a2-1.0) - np.power(L, -a2-1.0)) / a2
    S3 = 2.0 * u3 * (np.power(L, a3-1.0) - np.power(L, -a3-1.0)) / a3
    return S1+S2+S3


# --------------------------------------------------------------------------------
# Calculate the Ds
# --------------------------------------------------------------------------------
def compressibility(v, u1, a1, u2, a2, u3, a3):
    # This sum is what's in the ABQ manual (and what ABQ calculates with the data).
    # We get an error message which implies that u1 is what ABQ actually expects.
    # I believe the error message to be incorrect; setting u0=u1 typically results
    # in a much less compressible material, even though the error goes away.
    #  u0 = u1
    u0 = u1+u2+u3
    D1 = 3.0*(1.0-2.0*v) / (u0*(1.0+v))
    return [D1, 0.0, 0.0]
