# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 16:13:49 2023

@author: blr
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import copy

class xray:
    def __init__(self, n, Z, l, r, position=np.r_[0.0,0.0,0.0],direction=np.r_[0.0,0.0,1.0],energy=1.0):
        self.position=position
        self.direction=direction
        self.energy=energy
        self.n = n
        self.Z = Z
        self.terminated=False
        self.r = r
        self.l = l
        self.status="in cap"
    
    def location_relative_to_cap(self):
        if np.sum(self.position[0:2]**2) >= self.r**2 or self.position[2] < 0:
            self.status = "lost"
            return "lost"
        elif self.position[2] > self.l:
            self.status = "detected"
            return "detected"
        else:
            self.status = "in cap"
            return "in cap"
    
    def ef_compton(self,k,theta):
        return k / (1+k*(1-np.cos(theta)))
    
    def dsigma_dtheta_compton(self,k,theta):
        re = 2.8e-13
        ef = self.ef_compton(k, theta)
        return self.Z*re**2/2 * (ef/k)**2 * (ef/k + k/ef - np.sin(theta)**2)*2*np.pi*np.sin(theta)
    
    def sigma_compton(self,k):
        theta = np.linspace(0,np.pi,40)
        return np.trapz(self.dsigma_dtheta_compton(k,theta),theta)
    
    def sigma_photoelectric(self,k):
        re = 2.8e-13
        return 8*np.sqrt(2)*np.pi * re**2 * (1/137)**4 * self.Z**5 / k**3.5
    
    def distance_travelled(self,sigma):
        mfp = 1/self.n/sigma
        dist = scipy.stats.expon(scale=mfp).rvs()
        return dist
    
    def update_direction(self):
        thetas = np.linspace(0,np.pi,30)
        dsigma = self.dsigma_dtheta_compton(self.energy, thetas)
        theta = np.random.choice(thetas,p=dsigma/np.sum(dsigma))
        
        if self.direction[2]!=1.0:
            e2 = np.r_[self.direction[1],-self.direction[0],0.0]
        else:
            e2 = np.r_[1.0,0.0,0.0]
        e3 = np.cross(self.direction,e2)
        
        phi = np.random.rand(1)*2*np.pi
        dir_new = (np.cos(phi)*np.sin(theta)*e2 
                   +np.sin(phi)*np.sin(theta)*e3
                   +np.cos(theta)*self.direction)
        self.direction = dir_new/np.sqrt(np.sum(self.direction**2))
        self.energy = self.ef_compton(self.energy, theta)
        
        
        
    def update_position(self,x_travelled):
        self.position += x_travelled*self.direction
        
    
    def take_step(self):
        sigma_compton = self.sigma_compton(self.energy)
        sigma_photoelectric = self.sigma_photoelectric(self.energy)
        
        x_compton = self.distance_travelled(sigma_compton)
        x_photoelectric = self.distance_travelled(sigma_photoelectric)
        
        if x_compton > x_photoelectric:
            x_travelled = x_photoelectric
            self.status = "absorbed"
            self.terminated = True
        if x_compton < x_photoelectric:
            x_travelled = x_compton
            self.update_direction()
            
        self.update_position(x_travelled)
        where = self.location_relative_to_cap()
        if where in ["lost", "detected"]:
            self.terminated = True
        
        return x_compton, x_photoelectric, where
        
    def run_sim(self):
        energies = [copy.deepcopy(self.energy)]
        coords = [copy.deepcopy(self.position)]
        while self.terminated==False:
            self.take_step()
            coords.append(copy.deepcopy(self.position))
            energies.append(copy.deepcopy(self.energy))
            print(self.status)
            
        return coords, energies, self.status 
        
            
for i in range(2):
    x = xray(6.3e22,74.0,0.5,3.0,position=np.r_[0.0,0.0,0.1],direction=np.r_[0.0,0.0,1.0],energy=1.0)
    print(x.run_sim(),"\n\n")
        
        
        
        