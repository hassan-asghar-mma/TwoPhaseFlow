#!/usr/bin/env python3

# Purpose of this script is to investigate the volume ratios of spherical caps
# which share the same wetted area. This is the case for the growing droplet
# as long as it is within hysteresis range.
# NOTE: this volume ratio does not depend on the radius, but only on the
# two contact angles of the initial spherical cap and the one formed by the
# advancing contact angle.
from math import *
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300

# Initial contact angle of spherical cap in degree
theta0 = 60

# Angles for which to explicitly output the droplet volume
# STDOUT
expAngles = [80, 90, 110]

def volratio(theta0, thetaA):
    theta0 = theta0*pi/180.0    # Initial angle of spherical cap
    thetaA = thetaA*pi/180.0    # Advaning contact angle
    sinTerm = (sin(theta0)/sin(thetaA))**3
    cosNumerator = (2.0 + cos(thetaA))*(1.0 - cos(thetaA))**2
    cosDenominator = (2.0 + cos(theta0))*(1.0 - cos(theta0))**2

    return sinTerm*cosNumerator/cosDenominator

angles = [x for x in range(theta0, 121, 1)]
ratios = [volratio(theta0,x) for x in angles]

plt.plot(angles, ratios)
plt.xlabel("Contact angle in degree")
plt.ylabel("Multiples of initial volume")
plt.grid()
plt.savefig("contactAngle_dropletVolume.png")

for angle in expAngles:
    print("Volume ratio for initial angle", theta0, "and advancing angle", angle,
            "is", volratio(theta0, angle))
