#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

colors = ['#785EF0','#FE6100','#DC267F','#FFB000','#648FFF']


'''
k points analysis
a = 3.8034
'''
k_points = np.array([5, 7, 9, 11, 13, 15, 17, 19, 21])
k_energy = np.array([-7.12184310, -7.33284104, -7.26504444, -7.26333269, -7.27356074, -7.27243241, -7.27513999, -7.27245507, -7.27273996])
k_points = np.array([9, 11, 13, 15, 17, 19, 21])
k_energy = np.array([-7.26504444, -7.26333269, -7.27356074, -7.27243241, -7.27513999, -7.27245507, -7.27273996])
k_points = np.array([13, 15, 17, 19, 21])
k_energy = np.array([-7.27356074, -7.27243241, -7.27513999, -7.27245507, -7.27273996])
k = [x*1.60217663*10**(-22)*6.02214076*10**(23) for x in k_energy]

fig, ax = plt.subplots(1,1,figsize=(8,4))
ax.plot(k_points,k,linestyle='--',marker='o',color=colors[0])
ax.set_title('Finding the number of k points')
ax.set_xticks(k_points)
ax.set_xlabel('number of k points')
ax.set_ylabel('energy (kJ/mol)');