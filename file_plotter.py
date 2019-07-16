#file runner

import numpy as np
import matplotlib.pyplot as plt
import glob

file = []
datas = [];
us = [];
for f in glob.iglob("*.csv"): # generator, search immediate subdirectories   
    if f[0] == 'd':
        file.append(f)
        datas.append(np.loadtxt(f,  delimiter = ','))
    else:
        us.append(np.loadtxt(f,  delimiter = ','))


x = datas[0][0]
#initial condition
for j in range(len(datas)):
    f1 = plt.figure( figsize = (8,3) )
    plt.plot(x, datas[j][1])
    plt.plot(x, us[j][1])
    plt.legend(['Depth', 'Velocity'])
    plt.savefig('initial_'+file[j][10:15]+'.eps')
#final comparison of fluid depth
for j in range(len(datas)):
    f1 = plt.figure( figsize = (8,3) )
    plt.plot(x, datas[j][1])
    plt.plot(x, datas[j][-1])
    plt.legend(['t=0', 't=T'])
    plt.savefig('compare_'+file[j][10:15]+'.eps')
    
##f2 = plt.figure( figsize = (8,3) )
##plt.plot(x, datas[:][-1])
##plt.plot(x, us[:][-1])
##plt.legend(['Depth', 'Velocity'])
##plt.savefig('last.eps')
##f3 = plt.figure( figsize = (8,3) )
##plt.plot(x, datas[:][1])
##plt.plot(x, datas[:][-1])
##plt.legend(['t=0','t=20'])
##plt.savefig('compare.eps')

