from pylab import *
from scipy import integrate


def BB_lambda(wave,T):
    c = 3e8
    h = 6.626068e-34
    k = 1.3806503e-23

    temp = 2.*h*c**2./((wave*1e-6)**5)*1/(exp((h*c)/(wave*1e-6*k*T))-1)
    return temp

def BB_nu(nu,T):
    c = 3e8
    h = 6.626068e-34
    k = 1.3806503e-23

    temp = 2.*h*nu**3./(c**2)*1/(exp((h*nu)/(k*T))-1)
    return temp

ion()
clf()

color1=[0,0,1]
color2=[1.,0,1.]

wave_min_I=0
wave_max_I=122
#wave_min_I=120
#wave_max_I=240
itheta = 2
idusty = '6'
stepp_theta = 5.


flux,tau = loadtxt('slb_test_inten.s00'+idusty,usecols=[1,6],unpack=1)
tau = tau[wave_min_I:wave_max_I]
flux = flux[wave_min_I:wave_max_I]
wave,inten = loadtxt('slb_test_inten.i00'+idusty,usecols=[0,itheta],unpack=1)
wave = wave[wave_min_I:wave_max_I]
inten = inten[wave_min_I:wave_max_I]

sub1=subplot(211)
lwdt=5.
if (idusty=='1'):title(r'$T=1000,\tau$=0.01')
if (idusty=='2'):title(r'$T=1000,\tau$=0.1')
if (idusty=='3'):title(r'$T=1000,\tau$=1')
if (idusty=='4'):title(r'$T=1000,\tau$=10')
if (idusty=='5'):title(r'$T=1000,\tau$=100')
if (idusty=='6'):title(r'$T=1000,\tau$=1000')
plot(wave,(1-exp(-tau))*BB_lambda(wave,1000),'-',label='BB*(1-exp(-tau))',lw=lwdt,color=color1)
plot(wave,inten/(wave*1e-6),'--',label='dusty',lw=lwdt,color=color2)
#plot(wave,flux/wave*1e10)
ylabel('Intensity [W/m^2/sr/m]')
xlabel('wavelength [micron]')

integ = integrate.trapz(BB_lambda(wave,1000),wave*1e-6)
print integ
print  5.670373e-8*1e3**4/pi/integ

legend()
xscale('log')
yscale('log')
ylim(1e-1,1e10)
xlim(1e-1,1e4)


sub2=subplot(212)
plot(wave,1-(inten/(wave*1e-6))/((1-exp(-tau))*BB_lambda(wave,1000)),label='BB',lw=lwdt)
xlim(1e-1,1e4)
xscale('log')
ylim(-1.,1.)
ylabel('Difference (1-dusty/BB)')
xlabel('wavelength [micron]')

savefig('compareBB_'+idusty+'.eps')
