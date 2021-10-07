import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Solar House Model for QEA 3

# An an approximate model of the solar flux over time
def window_flux(t):
    flux = -361*math.cos(math.pi*t/(12*3600)) + 224*math.cos(math.pi*t/(6*3600)) + 210
    return flux # W / m^2

# Make "houses" that we can generate and optimize
class House:
    def __init__(self, height, width, depth, window, thickness_floor, thickness_wall):    
        # Shape of house
        self.height = height # m
        self.width = width # m
        self.depth = depth # m
        self.window = window # m
        self.thickness_floor = thickness_floor # m
        self.thickness_wall = thickness_wall # m

        self.construct_house()
        self.model_resistance()
        self.simulate_house()


    def construct_house(self):
        # Dependant sizes for proper window placement
        self.overhang = self.window/(np.tan(np.radians(75)) \
                        - np.tan(np.radians(25))) # m
        self.above_window = self.overhang*np.tan(np.radians(25)) # m
        self.below_window =  self.height-self.window-self.above_window # m

        # Floor Constant Properties
        DENSITY_FLOOR = 3000 #  kg/m^3
        SPECIFIC_HEAT_FLOOR = 800 # Joules / kg * K^-1

        # Floor Properties
        self.area_floor = self.width * self.depth * 2 # m^2
        self.capacity_floor = self.area_floor * self.thickness_floor \
                            * DENSITY_FLOOR * SPECIFIC_HEAT_FLOOR # m

        
        # Wall Parameters
        self.area_wall_inside = self.width*self.depth*2 + (self.above_window \
                                +self.below_window)*self.depth + self.height \
                                *self.depth + self.width*self.height*2 # m^2
        self.area_window = self.window * self.depth # m^2

        # Calculating the surface area of the outside of the house is a mess
        self.area_wall_outside = (self.width+2*self.thickness_wall) \
                                * (self.depth + 2*self.thickness_wall) \
                                + (self.height+2*self.thickness_wall)\
                                * (self.depth+2*self.thickness_wall) \
                                + ((self.width+self.overhang) \
                                + 2*self.thickness_wall) \
                                * (self.depth+2*self.thickness_wall) \
                                + ((self.above_window+self.below_window) \
                                + 2*self.thickness_wall) \
                                * (self.depth+2*self.thickness_wall) \
                                + (self.height+2*self.thickness_wall) \
                                + (self.width*2*self.thickness_wall) \
                                + self.overhang*self.thickness_wall*2
    
    def model_resistance(self):
        # Convection Coefficients
        H_INDOOR = 15 # W / m^2 * K^-1
        H_OUTDOOR = 30 # W / m^2 * K^-1
        H_WINDOW = 0.7 # W / m^2 * K^-1

        # Conduction coefficients
        K_FIBERGLASS = 0.04 # W / m * K^-1

        # Resistances
        r_1 = 1/(H_INDOOR * self.area_floor) # W/K
        r_2 = 1/(H_INDOOR * self.area_wall_inside) # W/K
        r_3 = 1/(H_INDOOR * self.area_window) # W/K
        r_4 = self.thickness_wall/(K_FIBERGLASS * self.area_wall_inside) # W/K
        r_5 = 1/(H_WINDOW * self.area_window) # W/K
        r_6 = 1/(H_OUTDOOR * self.area_wall_outside) # W/K
        r_7 = 1/(H_OUTDOOR * self.window) # W/K
        self.r_tot = r_1 + 1/(1/(r_2+r_4+r_6) + 1/(r_3+r_5+r_7)) # Total resistance W/K
        self.voltage_div = 1/(1/(r_2+r_4+r_6) + 1/(r_3+r_5+r_7))/self.r_tot

    def dTdt(self,t,T_floor):
        TEMP_OUT = -3 # C
        res = (window_flux(t)*self.area_window - (T_floor - TEMP_OUT)/self.r_tot)/self.capacity_floor # K/S
        return res

    def simulate_house(self):
        # Simulation Conditions
        TEMP_OUT = -3 # C
        T_IN_0 = 21 # C
        DAYS = 100 # d
        T_0_EVAL = 80 # time to stabilize d
        T_SPAN = [0,86400*DAYS] # s
        T_SPAN_EVAL = [86400*T_0_EVAL,86400*DAYS] # s
        self.results = solve_ivp(self.dTdt, T_SPAN, [T_IN_0], t_eval = range(86400*DAYS))
        self.air_temp = (self.results.y.T - TEMP_OUT) * self.voltage_div + TEMP_OUT
        self.std = np.std(self.air_temp[T_SPAN_EVAL[0]:T_SPAN_EVAL[1]])
        self.avg = np.mean(self.air_temp[T_SPAN_EVAL[0]:T_SPAN_EVAL[1]])

    def graph_house(self):
        plt.plot(self.results.t/86400, self.air_temp)
        plt.xlabel("Time (days)")
        plt.ylabel("Temp (C)")
        plt.title("Air Temp over time")
        plt.show()


test_house = House(3,5.1,5,2.6,.3,.031)
test_house.graph_house()
print([test_house.std, test_house.avg])