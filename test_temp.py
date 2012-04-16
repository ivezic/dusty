from pylab import *

ion()
clf()

#y_old,T_old = loadtxt('/home/fheymann/work/dusty_benchmark/ivezic97/ivezic97.r004',unpack=1,usecols=[0,1],skiprows=4)
y_new,fde,fds,fsbol,T_new = loadtxt('fort.7777',unpack=1,usecols=[0,1,2,3,4])
last = y_new.size-1
max_n = last
print y_new[0],y_new[1],y_new.size
y_new = arange(y_new.size)*100./y_new.size
plot(y_new[:max_n],fde[:max_n],'|-',label='fde')
plot(y_new[:max_n],fds[:max_n],'|-',label='fds')
plot(y_new[:max_n],fsbol[:max_n],'|-',label='fsbol')
plot(y_new[:max_n],fde[:max_n]+fds[:max_n]+fsbol[:max_n],'x-')

fbol = fde+fds+fsbol

print(fde[last]+fds[last]+fsbol[last])
legend()
#print max(fde+fds+fsbol)/min(fde+fds+fsbol)
#clf()
#plot(y_new,fde)
#n=1000
#y_new = arange(n)*1./(n-1)*1000.
#plot(y_new,4.+1./((y_new/100)**2))

#xlim(0,100)
#ylim(8,13)
#ylim(0,100)
#xlim(1,2)
#ylim(10,14)
#xscale('log')
#yscale('log')

#clf()
#plot(y_new,T_new,'x')
#xscale('log')

#clf()
#plot(y_new[:y_new.size-2],fbol[:y_new.size-2]/fbol[1:y_new.size-1],'x')
#xscale('log')
