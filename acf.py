#!/usr/bin/env python

import sys
from numpy import *
from pylab import *

T  = 40.0
sf = 300.0
t = arange(0,T,1.0/sf)

n  = int(T*sf)

f=zeros(n)
g=zeros(n)

# Time-domain function
freq=0.05
w = 2*pi*freq
#f = (0.0 + sin(w*t)) * sin(10*w*t)
#f = 10.0 * sin(10*w*t) + 1.0 * sin(w*t)
f = 5.0*sin(2*pi*0.5*t)

# Fourier transform
s=arange(-n/2,n/2)/T
tmp = fft(f)/n
g[:n/2] = abs(tmp[n/2:])
g[n/2:] = abs(tmp[:n/2])

fig1 = figure(1,figsize=(8,12))
subplot(4,1,1)
plot(t,f)

subplot(4,1,2)
thd=0.01
for i in range(n):
    if g[i]>thd:
        plot([s[i],s[i]],[0,g[i]],linewidth=2.0,color='k')
setp(gca(),xlim=[-1,1])

subplot(4,1,3)
a=correlate(f,f,mode='same')
plot(range(-n/2,n/2),a)
setp(gca(),xlim=[-1000,1000],yticks=[])

subplot(4,1,4)
b=zeros(n)
tmp = fft(a)
tmax=max(tmp)
tmp=tmp/tmax
b[:n/2] = abs(tmp[n/2:])
b[n/2:] = abs(tmp[:n/2])
for i in range(n):
    if b[i]>thd:
        plot([s[i],s[i]],[0,b[i]],linewidth=2.0,color='k')
setp(gca(),xlim=[-1,1])

savefig('foo.pdf')
