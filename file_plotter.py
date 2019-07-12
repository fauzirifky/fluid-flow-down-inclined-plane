
#file runner

import numpy as np
import matplotlib.pyplot as plt

datas = np.loadtxt('data_plot.csv',  delimiter = ',')
us = np.loadtxt('u_plot.csv',  delimiter = ',')
x = datas[:][0]
f1 = plt.figure( figsize = (8,3) )
plt.plot(x, datas[:][1])
plt.plot(x, us[:][1])
plt.legend(['Depth', 'Velocity'])
plt.savefig('initial.eps')
f2 = plt.figure( figsize = (8,3) )
plt.plot(x, datas[:][-1])
plt.plot(x, us[:][-1])
plt.legend(['Depth', 'Velocity'])
plt.savefig('last.eps')
f3 = plt.figure( figsize = (8,3) )
plt.plot(x, datas[:][1])
plt.plot(x, datas[:][-1])
plt.legend(['t=0','t=20'])
plt.savefig('compare.eps')

