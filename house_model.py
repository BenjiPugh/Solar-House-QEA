import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Solar House Model for QEA 3

# Make "houses" that we can generate and optimize
class House:
    def __init__(self, height, width, depth, window):    
        # Shape of house
        self.height = height # m
        self.width = width # m
        self.depth= depth # m
        self.window = window # m
        self.construct_house()

    def construct_house(self):
        self.overhang = self.window/(np.tan(np.adians(75)) \
                        - np.tan(np.radians(25))) # m
        self.above_window = self.overhang*np.tan(np.radians(25)) # m
        self.below_window =  self.height-self.window-self.above_window # m


# Floor Parameters

area_floor = 5.1 * 5 * 2 # m^2
density_floor = 3000 #  kg/m^3
thickness_floor = 0.1 # meters
specific_heat_floor = 800 # Joules / kg * K^-1
capacity_floor = area_floor * thickness_floor * density_floor \
                * specific_heat_floor

# Convection Coefficients
h_indoor = 15 # W / m^2 * K^-1
h_outdoor = 30 #    W / m^2 * K^-1

# Wall Parameters
wall_thickness = 0.1 # meters
area_wall_inside = width*depth*2 + (above_window+below_window)*depth + height*depth + width*height*2 # m^2
area_window = window * depth # m^2
area_wall_outside = (width+2*wall_thickness) * (depth + 2*wall_thickness) \
                    + (height+2*wall_thickness)*(depth+2*wall_thickness) \
                    + ((width+overhang)+2*wall_thickness)*(depth+2*wall_thickness) \
                    + ((above_window+below_window)+2*wall_thickness)*(depth+2*wall_thickness) \
                    + (height+2*wall_thickness)+(width*2*wall_thickness) \
                    + overhang*wall_thickness*2

# Conduction coefficients
k_fiberglass = 0.04 # W / m * K^-1
h_window = 0.7 # W / m^2 * K^-1

# Resistances
r_1 = 1/(h_indoor * area_floor) # W/K
r_2 = 2/(h_indoor * (area_wall_inside + area_window)) # W/K
r_3 = wall_thickness/(k_fiberglass * area_wall_inside) # W/K
r_4 = 1/(h_window * area_window) # W/K
r_5 = 1/(h_outdoor * (area_wall_outside + area_window)) # W/K
r_tot = r_1 + r_2 + (1/(1/r_3+1/r_4)) + r_5 # Total resistance W/K

# Initial Temperatures
TEMP_OUT = -3 # C
T_IN_0 = 21 # C


# An an approximate model of the solar flux over time
def window_flux(t):
    flux = -361*math.cos(math.pi*t/(12*3600)) + 224*math.cos(math.pi*t/(6*3600)) + 210
    return flux # W / m^2

def dTdt(t,T_floor):
    res = (window_flux(t)*area_window - (T_floor - TEMP_OUT)/r_tot)/capacity_floor # K/S
    return res

voltage_div = (r_2 + (1/(1/r_3+1/r_4)) + r_5) / (r_tot)


# Find numerical solution to ODE
days = 40
t_span = [0,86400*days] # s
results = solve_ivp(dTdt, t_span, [T_IN_0], t_eval = range(86400*days))
print(results)
print(type(results))

air_temp = (results.y.T - TEMP_OUT) * voltage_div + TEMP_OUT
plt.plot(results.t/86400, air_temp)
plt.xlabel("Time (days)")
plt.ylabel("Temp (C)")
plt.title("Air Temp over time")
plt.show()