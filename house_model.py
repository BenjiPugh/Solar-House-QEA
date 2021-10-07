import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Solar House Model for QEA 3

# An an approximate model of the solar flux over time
def window_flux(t):
    flux = -361*math.cos(math.pi*t/(12*3600)) + 224*math.cos(math.pi*t/(6*3600)) + 210
    return flux # W / m^2

# Functions to plot min and max
def annot_max(x,y, ax=None):
    xmax = x[np.argmax(y)]
    ymax = y.max()
    text= "Maximum Temperature\nx={:.3f}, y={:.3f}".format(xmax, ymax)
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="axes fraction",
              arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(xmax, ymax), xytext=(0.94,0.8), **kw)

def annot_min(x,y, ax=None):
    xmin = x[np.argmin(y)]
    ymin = y.min()
    text= "Minimum Temperature\nx={:.3f}, y={:.3f}".format(xmin, ymin)
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="axes fraction",
              arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(xmin, ymin), xytext=(0.94,0.2), **kw)


# Simulation Conditions
SEC_DAY = 86400 # s
TEMP_OUT = -3 # C
T_IN_0 = 21 # C
DAYS = 100 # d
T_0_EVAL = 80 # time to stabilize d
T_SPAN = [0,SEC_DAY*DAYS] # s
T_SPAN_EVAL = [SEC_DAY*T_0_EVAL,SEC_DAY*DAYS] # s

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
        res = (window_flux(t)*self.area_window - (T_floor - TEMP_OUT)/self.r_tot) \
            / self.capacity_floor # K/S
        return res

    def simulate_house(self):
        self.results = solve_ivp(self.dTdt, T_SPAN, [T_IN_0], t_eval = range(SEC_DAY*DAYS))
        self.air_temp = (self.results.y.T - TEMP_OUT) * self.voltage_div + TEMP_OUT
        self.std = np.std(self.air_temp[T_SPAN_EVAL[0]:T_SPAN_EVAL[1]])
        self.avg = np.mean(self.air_temp[T_SPAN_EVAL[0]:T_SPAN_EVAL[1]])

    def graph_simulation(self):
        plt.plot(self.results.t/SEC_DAY, self.air_temp)
        plt.xlabel("Time (days)")
        plt.ylabel("Temp (C)")
        plt.title("Air Temperature Over Time")
        plt.show()

    def graph_eval(self):
        plt.plot(self.results.t[T_SPAN_EVAL[0]:T_SPAN_EVAL[1]]/SEC_DAY \
               , self.air_temp[T_SPAN_EVAL[0]:T_SPAN_EVAL[1]])
        plt.xlabel("Time (s)")
        plt.ylabel("Temp (C)")
        plt.title("Air Temperature Over Evaluation Period")
        plt.show()

    def eval_day(self):
        self.day_time = (self.results.t[T_0_EVAL*SEC_DAY:T_0_EVAL*SEC_DAY+SEC_DAY] \
                        /SEC_DAY*24-0)-T_0_EVAL*24
        self.day_air_temp = self.air_temp[T_0_EVAL*SEC_DAY:T_0_EVAL*SEC_DAY+SEC_DAY]
        self.day_max = np.max(self.day_air_temp)
        self.day_min = np.min(self.day_air_temp)
        self.day_ave = np.mean(self.day_air_temp)
        self.day_sd = np.std(self.day_air_temp)

    def graph_day(self):
        self.eval_day()
        print(self.day_time)
        plt.plot(self.day_time, self.day_air_temp)
        # day_sd_vis = [self.day_ave+self.day_sd,self.day_ave-self.day_sd] not necessary
        annot_max(self.day_time, self.day_air_temp)
        annot_min(self.day_time, self.day_air_temp)
        plt.axhline(self.day_ave,color="r",linestyle="--")
        kw = dict(xycoords='data',textcoords="axes fraction",ha="left", va="top")
        plt.annotate("Average Temperature", xy=(12, self.day_ave), xytext=(0.02,0.55), **kw)
        plt.xlabel("Time (Hours)")
        plt.ylabel("Temp (C)")
        plt.title("Air Temperature For One Day")
        plt.show()


test_house = House(3,5.1,5,2.6,.3,.031)
test_house.graph_day()
print([test_house.std, test_house.avg])