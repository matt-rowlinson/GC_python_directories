#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 15:01:45 2018

@author: ee11mr
"""
import numpy as np

##initialise variables 
n_lons=720
n_lats=360
tc_lons=np.arange(0.25,360.25,.5)
tc_lats=np.arange(-89.25,90.25,0.5)
lons_edges=np.arange(0,360.5,.5)
lats_edges=np.arange(-90.0,90.5,.5)

R=6371000.0    #radius of earth in km
mdi=-999.99
surface_area_earth=np.zeros((n_lons,n_lats))+mdi
#lon lats in radians
lons_rad=lons_edges*(2*np.pi/360.0)
lats_rad=lats_edges*(2*np.pi/360.0)

#calculate global surface area on lon lat grid
#surface area = R^2(lambda1-lambda2)(sin(phi1)-sin(phi2)) - lambda in radians, phi = -pi/2 to pi/2
for ilon in range(n_lons):
    for ilat in range(n_lats):
        #get terms
        term1=R**2
        term2=lons_rad[ilon+1]-lons_rad[ilon]
        term3=np.sin(lats_rad[ilat+1])-np.sin(lats_rad[ilat])
        #;surface area
        tmp_sa=term1*term2*term3
        surface_area_earth[ilon,ilat]=tmp_sa
#should be roughly 510 million km2
#surface_area_earth=np.reshape(surface_area_earth,64800)
print('Surface Area of Earth (km^2)= ',np.sum(surface_area_earth))
