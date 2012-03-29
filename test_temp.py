from pylab import *

ion()
#clf()

#y_old,T_old = loadtxt('/home/fheymann/work/dusty_benchmark/ivezic97/ivezic97.r004',unpack=1,usecols=[0,1],skiprows=4)
y_new,fde,fds,fsbol,T_new = loadtxt('fort.7777',unpack=1,usecols=[0,1,2,3,4])
plot(y_new,fde)
plot(y_new,fds)
plot(y_new,fsbol)
plot(y_new,fde+fds+fsbol)
print max(fde+fds+fsbol)/min(fde+fds+fsbol)
