"""Angular map processing for EasyPizza, kept free of Tkinter so it can be tested headless.

Conventions (unchanged from the original program):
- theta is the polar emission angle in degrees; it is the radius of the map.
- phi is the azimuth: read in degrees, returned in radians, shifted so the map starts at
  zero and cropped at 2*pi.
- .xpd files hold "Intensity Theta Phi"; each theta ring is normalised as
  chi = (I - <I>) / <I>. .hol files use columns 0, 1, 4 and .dat files columns 0, 1, 2
  as (theta, phi, chi).
- Symmetry modes: 0 raw data, 1 60 degrees, 2 120 degrees, 3 90 degrees.
"""
import numpy as np
from scipy import interpolate

RAW, SYM60, SYM120, SYM90 = 0, 1, 2, 3


def load_map(path, symmetry=RAW, rotate=0):
    """Read a map file and return the scattered points (THETA, PHI, ki) that get interpolated."""
    ext = path.split('.')[-1]
    A = np.loadtxt(path, skiprows=1)

    if ext == "xpd":
        THETA, PHI, ki = _xpd_points(A[:, 0], A[:, 1], symmetry)
    elif ext == "hol":
        THETA, PHI, ki = A[:, 0], A[:, 1], A[:, 4]
    elif ext == "dat":
        THETA, PHI, ki = A[:, 0], A[:, 1], A[:, 2]
    else:
        raise ValueError("Unsupported file type .%s (expected .xpd, .hol or .dat)" % ext)

    PHI = np.array(PHI)
    #Shift to get the XPD pizza centralized
    PHI = np.pi*PHI/180 - np.pi*rotate/180

    PHI_, THETA_, ki_ = PHI, THETA, ki
    copies = {SYM60: (6, 60), SYM120: (3, 120), SYM90: (4, 90)}
    if symmetry in copies:
        count, step = copies[symmetry]
        for i in range(1, count):
            PHI = np.append(PHI, PHI_ + i*step*np.pi/180)
            THETA = np.append(THETA, THETA_)
            ki = np.append(ki, ki_)

    PHI = PHI - min(PHI)  #Shifts PHI to zero
    del_index = np.where(PHI > 2*np.pi)  #Index of everything greater than 2pi
    THETA = np.delete(THETA, del_index)  #Deletes indexes
    ki = np.delete(ki, del_index)
    PHI = np.delete(PHI, del_index)
    return THETA, PHI, ki


def _xpd_points(intensity, theta_degree, symmetry):
    """Normalise each theta ring and, for the symmetry modes, fold it into one sector."""
    ki = np.array([])
    THETA, PHI = [], []

    for t in range(int(min(theta_degree)), int(max(theta_degree))+1, 3):
        I_t = intensity[np.where(theta_degree == t)]
        I_tm = []
        if symmetry == RAW:
            for i in range(len(I_t)):
                THETA.append(t)
                PHI.append(3*i)
            I_m = np.sum(I_t)/len(I_t)
            NI = (I_t-I_m)/I_m
            ki = np.append(ki, NI)

        if symmetry == SYM60:
            if len(I_t) % 4 == 0:
                for i in range(len(I_t)//4):
                    I_tm.append((I_t[i]+I_t[i+len(I_t)//4]+I_t[i+2*len(I_t)//4]+I_t[i+3*len(I_t)//4])/2.)
                    THETA.append(t)
                    PHI.append(3*i)
            else:
                for i in range(len(I_t)//4+1):
                    if I_t[i] == I_t[-1]:
                        I_tm.append(I_t[i])
                        THETA.append(t)
                        PHI.append(3*i)
                    else:
                        I_tm.append((I_t[i]+I_t[i+len(I_t)//4]+I_t[i+2*len(I_t)//4]+I_t[i+3*len(I_t)//4])/2.)
                        THETA.append(t)
                        PHI.append(3*i)
            I_m = np.sum(I_tm)/len(I_tm)
            NI = (I_tm - I_m)/I_m
            ki = np.append(ki, NI)

        if symmetry == SYM120:
            if len(I_t) > 230/3.:
                if len(I_t) % 2 == 0:
                    for i in range(len(I_t)//2):
                        I_tm.append((I_t[i]+I_t[i+len(I_t)//2])/2.)
                        THETA.append(t)
                        PHI.append(3*i)
                else:
                    for i in range(len(I_t)//2+1):
                        if I_t[i] == I_t[-1]:
                            I_tm.append(I_t[i])
                            THETA.append(t)
                            PHI.append(3*i)
                        else:
                            I_tm.append((I_t[i]+I_t[i+len(I_t)//2])/2.)
                            THETA.append(t)
                            PHI.append(3*i)
            else:
                for i in range(len(I_t)):
                    I_tm.append(I_t[i])
                    THETA.append(t)
                    PHI.append(3*i)
            I_m = np.sum(I_tm)/len(I_tm)
            NI = (I_tm - I_m)/I_m
            ki = np.append(ki, NI)

        if symmetry == SYM90:
            if len(I_t) % 3 == 0:
                for i in range(len(I_t)//3):
                    I_tm.append((I_t[i]+I_t[i+len(I_t)//3]+I_t[i+2*len(I_t)//3])/2.)
                    THETA.append(t)
                    PHI.append(3*i)
            else:
                for i in range(len(I_t)//3+1):
                    if I_t[i] == I_t[-1]:
                        I_tm.append(I_t[i])
                        THETA.append(t)
                        PHI.append(3*i)
                    else:
                        I_tm.append((I_t[i]+I_t[i+len(I_t)//3]+I_t[i+2*len(I_t)//3])/2.)
                        THETA.append(t)
                        PHI.append(3*i)
            I_m = np.sum(I_tm)/len(I_tm)
            NI = (I_tm - I_m)/I_m
            ki = np.append(ki, NI)

    return THETA, PHI, ki


def interpolate_map(THETA, PHI, ki, symmetry=RAW):
    """Interpolate the scattered points on a 600 x 600 grid; returns (gtheta, gphi, z)."""
    if symmetry == RAW:
        gx, gy = np.mgrid[min(THETA):max(THETA):600j, min(PHI):max(PHI):600j]
        z = interpolate.griddata((THETA, PHI), ki, (gx, gy), method='cubic', fill_value=ki.min())
        return gx, gy, z

    xx = THETA*np.cos(PHI)
    yy = THETA*np.sin(PHI)

    gx, gy = np.mgrid[min(xx):max(xx):600j, min(yy):max(yy):600j]

    z = interpolate.griddata((xx, yy), ki, (gx, gy), method='cubic', fill_value=ki.min())
    gphi = np.arctan2(gy, gx)
    gtheta = np.sqrt(gx**2 + gy**2)

    gphi = gphi - np.amin(gphi)

    h, l = np.where(gtheta < min(THETA))
    z[h, l] = min(ki)
    return gtheta, gphi, z
