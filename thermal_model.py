import numpy as np
import matplotlib.pyplot as plt
from model import*

class thermal_model:
    """Class for creating thermal models."""
    def __init__(self, name, depth):
        """Name is name of the model and depth is starting depth of the moho."""
        self.name=name
        self.original_z=np.linspace(0, depth, 100)
        self.time=0
        self.depth=depth
        self.models=[]
    
    def calculate_geotherm(self, T_m, A, L, k, v):
        """Function for calculating initial geotherm. T_m is temperature of moho, A is relative heat production,
        L is moho depth, k is average conductivity above moho and v is denudation(positive) or errosion(negative)."""
        geotherm=[]
        Pe=v*L/k
        for z in self.original_z:
            T=-2*A*z/(Pe*L)+(1+2*A/Pe)*((1-np.exp(-Pe*z/L))/(1-np.exp(-Pe)))
            geotherm.append(T*T_m)
        self.geotherm=geotherm
        self.T_m=T_m
        self.time=0
        self.models.append(model(self.original_z, self.geotherm, self.time))
        self.A=A
        self.k=k
        self.L=L
    
    def collide(self, v, thickening, time, step=None):
        """Function for calculating futher geotherms. v is denudation(positive) or errosion(negative),
        thickening is absolute thickening of the model in time and step is amount of steps, if None only 1 step is used"""
        if step==None:
            step=time
        aim_time=self.time+time
        depth=self.depth
        while(self.time<aim_time):
            self.time+=step
            depth+=thickening/time*step-v*step
            effect=v-thickening/time
            z_list=np.linspace(0, depth, 100)
            geotherm=[]
            Pe=effect*depth/self.k
            for z in z_list:
                T=-2*self.A*z/(Pe*depth)+(1+2*self.A/Pe)*((1-np.exp(-Pe*z/depth))/(1-np.exp(-Pe)))
                geotherm.append(T*self.T_m)
            self.models.append(model(z_list, geotherm, self.time))
        self.depth=depth
    
    def plot_models(self, ylim=100):
        plt.figure(figsize=(10, 10))
        """Helpper function for plotting calculated models."""
        for i in self.models:
            plt.plot(i.geotherm, i.z, label=i.time)
            print(i.z[-1])
        plt.ylim([ylim, 0])
        plt.legend()
        plt.xlim([0, self.T_m+50])
        plt.grid()
        plt.ylabel('Depth(km)')
        plt.xlabel('Temperature (°C)')