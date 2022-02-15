import numpy as np
import random
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit

def configPlot_time(s,L,f,i,n_):
	""" This module plots the configuration """
	X, Y = np.meshgrid(range(L), range(L))
	config = s.reshape((L, L))
	sp =  f.add_subplot(3, 3, n_ )
	plt.setp(sp.get_yticklabels(), visible=False)
	plt.setp(sp.get_xticklabels(), visible=False)      
	plt.pcolormesh(X, Y, config)
	plt.title('Time=%d'%i); plt.axis('tight')
	return

def func(x, a):
	return np.exp(-x/a)

def main():
	L=50
	N=L*L
	nmcstep=10000
	
	""" energy, manetization """
	f = open('output3.dat','r')
	lines = f.readlines()
	f.close()
	x = []
	e = []
	m = []
	#ma= []
	for line in lines:
		p=line.split()
		x.append(int(p[0]))
		e.append(float(p[1])/N)
		m.append(float(p[2]))
		#ma.append(float(p[3]))
	fig=plt.figure(1,figsize=(6,5),dpi=100,facecolor='w',edgecolor='w')
	ax1 = fig.add_subplot(111)
	ax1.plot(x,e,'r-',markersize=3,label='energy density')
	ax1.plot(x,m,'b-',markersize=3,label='magnetization')
	#ax1.plot(x,ma,'k-',markersize=3,label='magnetization exact')
	ax1.set_xlabel(r'$MC \, time \, steps$')
	ax1.legend()
	plt.show()
	
	""" trajectory snapshots """
	f = open('output4.dat','r')
	lines = f.readlines()
	f.close()
	s=[]
	s1=[]
	s2=[]
	s3=[]
	s4=[]
	s5=[]
	for line in lines:
		p=line.split()
		s.append(float(p[0]))
		s1.append(float(p[1]))
		s2.append(float(p[2]))
		s3.append(float(p[3]))
		s4.append(float(p[4]))
		s5.append(float(p[5]))
	s=np.asarray(s)
	s1=np.asarray(s1)
	s2=np.asarray(s2)
	s3=np.asarray(s3)
	s4=np.asarray(s4)
	s5=np.asarray(s5)
	N=len(s)
	L=int(np.sqrt(N))

	f = plt.figure(figsize=(15, 15), dpi=80)

	configPlot_time(s, L, f, 0, 1)
	configPlot_time(s1, L, f, 1*nmcstep, 2)
	configPlot_time(s2, L, f, 5*nmcstep, 3)
	configPlot_time(s3, L, f, 15*nmcstep, 4)
	configPlot_time(s4, L, f, 50*nmcstep, 5)
	configPlot_time(s5, L, f, 100*nmcstep, 6)

	plt.show()
	
	
	""" Magnetization autocorrelation function """
# 	f = open('output5.dat','r')
# 	lines = f.readlines()
# 	f.close()
# 	x = []
# 	y = []
# 	for line in lines[0:200]:
# 		p=line.split()
# 		x.append(int(p[0]))
# 		y.append(float(p[1]))
# 	popt, pcov = curve_fit(func, x, y, p0=1000.)
# 	ys=[]
# 	for i in range(len(x)):
# 		ys.append(func(x[i],popt[0]))
# 	
# 	fig=plt.figure(figsize=(6,5),dpi=100,facecolor='w',edgecolor='w')
# 	ax1 = fig.add_subplot(111)
# 	ax1.plot(x,y,'ro-',markersize=3)
# 	ax1.plot(x,ys,'k--',markersize=3,label=r'$\exp(-t/\tau)$')
# 	ax1.set_xlabel(r'$t \; [MC \, time \, steps]$')
# 	ax1.set_ylabel(r'$\chi(t)$')
# 	ax1.legend()
# 	ax1.text(0.5,0.82,(r'$\tau =$ %4.1f $\pm$ %4.1f' % (popt[0],pcov[0][0])),fontsize=10,transform=ax1.transAxes)
# 	plt.show()
# 	
# 	return
	
main()

