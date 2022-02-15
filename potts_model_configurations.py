# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 10:04:10 2022

@author: louis
"""

import numpy as np
import random
import matplotlib.pyplot as plt


def configPlot(s,L):
	""" This module plots the configuration """
	f = plt.figure(figsize=(5, 5), dpi=100)
	config = s.reshape((L, L))
	sp =  f.add_subplot(111)
	X, Y = np.meshgrid(range(L), range(L))
	plt.setp(sp.get_yticklabels(), visible=False)
	plt.setp(sp.get_xticklabels(), visible=False)      
	plt.pcolormesh(X, Y, config)
	plt.show()
	return

def main():
	f = open('output1.dat','r')
	lines = f.readlines()
	f.close()
	s=[]
	for line in lines:
		p=line.split()
		s.append(float(p[0]))
	s=np.asarray(s)
	N=len(s)
	L=int(np.sqrt(N))
	configPlot(s,L)
	return
	
main()