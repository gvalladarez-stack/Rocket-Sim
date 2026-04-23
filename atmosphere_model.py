'''Compute ICAO Standard Atmosphere properties at a single altitude.
Returns temperature, pressure, density, and speed of sound.'''
from math import *
### Assumptions ###
# Conversion between km and ft is precisely 3280.84
# Air is treated as an ideal gas
# The atmosphere is modeled using the ICAO Standard Atmosphere
# Acceleration due to gravity is considered constant with varying altitude
# Earth's radius is constant
# Valid for geometric altitudes from 0 to 105 km inclusive.
# The function converts geometric altitude to geopotential altitude internally.
# Effects of wind, humidity, solar activity, etc are neglected

# Given constants (in English units)
P_sl = 2116.2 
p_sl = .002377
t_sl = 518.69
r_earth = 20925650
g = 32.2
R_air = 1716
y = 1.4

# Where Atmosphere Changes Type 
tropo_s = 36089
tropo_p = 82349
strato_s = 47*3280.84
strato_p = 53*3280.84
meso_s = 79*3280.84
meso_p = 90*3280.84
thermo_s = 105*3280.84

MAX_ALT_FT = thermo_s

def atmosphere(alt_ft):
    if alt_ft < 0 or alt_ft > MAX_ALT_FT:
        raise ValueError(
            f"Altitude must be between 0 and {MAX_ALT_FT:.0f} ft "
            f"(model valid to 105 km geometric altitude)"
        )
    alt_potential = r_earth * alt_ft / (r_earth + alt_ft)

    if alt_potential <= tropo_s:  # Troposphere
        a = -.0065 * 1.8 / 3.28084
        temp = a * alt_potential + t_sl
        Press = P_sl * (temp / t_sl) ** (-g / (a * R_air))
        rho = p_sl * (temp / t_sl) ** (-((g / (a * R_air)) + 1))

    elif alt_potential > tropo_s:
        a = -.0065 * 1.8 / 3.28084
        temp = 216.66 * 1.8
        Press = P_sl * (temp / t_sl) ** (-g / (a * R_air))
        rho = p_sl * (temp / t_sl) ** (-((g / (a * R_air)) + 1))

        if alt_potential <= tropo_p:
            Press = Press * e ** (-g * (alt_potential - tropo_s) / (R_air * temp))
            rho = rho * e ** (-g * (alt_potential - tropo_s) / (R_air * temp))

        else:
            Press = Press * e ** (-g * (tropo_p - tropo_s) / (R_air * temp))
            rho = rho * e ** (-g * (tropo_p - tropo_s) / (R_air * temp))

        if alt_potential > tropo_p and alt_potential <= strato_s:
            a = .003 * 1.8 / 3.28084
            temp = (alt_potential - tropo_p) * a + temp
            Press = Press * (temp / (216.66 * 1.8)) ** (-g / (a * R_air))
            rho = rho * (temp / (216.66 * 1.8)) ** (-((g / (a * R_air)) + 1))

        if alt_potential > strato_s and alt_potential <= strato_p:
            temp = 282.66 * 1.8
            Press = Press * (temp / (216.66 * 1.8)) ** (-g / ((.003 * 1.8 / 3.28084) * R_air))
            Press = Press * e ** (-g * (alt_potential - strato_s) / (R_air * temp))
            rho = rho * (temp / (216.66 * 1.8)) ** (-((g / ((.003 * 1.8 / 3.28084) * R_air)) + 1))
            rho = rho * e ** (-g * (alt_potential - strato_s) / (R_air * temp))

        if alt_potential > strato_p and alt_potential <= meso_s:
            a = -.0045 * 1.8 / 3.28084
            temp = 282.66 * 1.8
            temp = (alt_potential - strato_p) * a + temp
            Press = Press * (282.66 * 1.8 / (216.66 * 1.8)) ** (-g / ((.003 * 1.8 / 3.28084) * R_air))
            Press = Press * e ** (-g * (strato_p - strato_s) / (R_air * 282.66 * 1.8))
            Press = Press * (temp / (282.66 * 1.8)) ** (-g / (a * R_air))
            rho = rho * (282.66 * 1.8 / (216.66 * 1.8)) ** (-((g / ((.003 * 1.8 / 3.28084) * R_air)) + 1))
            rho = rho * e ** (-g * (strato_p - strato_s) / (R_air * 282.66 * 1.8))
            rho = rho * (temp / (282.66 * 1.8)) ** (-((g / (a * R_air)) + 1))

        if alt_potential > meso_s and alt_potential <= meso_p:
            temp = 165.66 * 1.8
            Press = Press * (282.66 * 1.8 / (216.66 * 1.8)) ** (-g / ((.003 * 1.8 / 3.28084) * R_air))
            Press = Press * e ** (-g * (strato_p - strato_s) / (R_air * 282.66 * 1.8))
            Press = Press * (165.66 * 1.8 / (282.66 * 1.8)) ** (-g / ((-.0045 * 1.8 / 3.28084) * R_air))
            Press = Press * e ** (-g * (alt_potential - meso_s) / (R_air * 165.66 * 1.8))
            rho = rho * (282.66 * 1.8 / (216.66 * 1.8)) ** (-((g / ((.003 * 1.8 / 3.28084) * R_air)) + 1))
            rho = rho * e ** (-g * (strato_p - strato_s) / (R_air * 282.66 * 1.8))
            rho = rho * (temp / (282.66 * 1.8)) ** (-((g / ((-.0045 * 1.8 / 3.28084) * R_air)) + 1))
            rho = rho * e ** (-g * (alt_potential - meso_s) / (R_air * temp))

        if alt_potential > meso_p and alt_potential <= thermo_s:
            a = .004 * 1.8 / 3.28084
            temp = 165.66 * 1.8
            temp = (alt_potential - meso_p) * a + temp
            Press = Press * (282.66 * 1.8 / (216.66 * 1.8)) ** (-g / ((.003 * 1.8 / 3.28084) * R_air))
            Press = Press * e ** (-g * (strato_p - strato_s) / (R_air * 282.66 * 1.8))
            Press = Press * (165.66 * 1.8 / (282.66 * 1.8)) ** (-g / ((-.0045 * 1.8 / 3.28084) * R_air))
            Press = Press * e ** (-g * (meso_p - meso_s) / (R_air * 165.66 * 1.8))
            Press = Press * (temp / (165.66 * 1.8)) ** (-g / (a * R_air))
            rho = rho * (282.66 * 1.8 / (216.66 * 1.8)) ** (-((g / ((.003 * 1.8 / 3.28084) * R_air)) + 1))
            rho = rho * e ** (-g * (strato_p - strato_s) / (R_air * 282.66 * 1.8))
            rho = rho * (165.66 * 1.8 / (282.66 * 1.8)) ** (-((g / ((-.0045 * 1.8 / 3.28084) * R_air)) + 1))
            rho = rho * e ** (-g * (meso_p - meso_s) / (R_air * 165.66 * 1.8))
            rho = rho * (temp / (165.66 * 1.8)) ** (-((g / (a * R_air)) + 1))

    c = (R_air * temp * y) ** 0.5
    # English to SI
    temp_SI = temp * (5/9)                 # Kelvin
    pressure_SI = Press * 47.8803          # Pascals
    rho_SI = rho * 515.379                 # kg/m³
    c_SI = c * 0.3048                      # m/s

    return temp_SI, pressure_SI, rho_SI, c_SI   
