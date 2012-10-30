#!/usr/bin/env python

import sys
from numpy import *
from pylab import *

C = array(map(float,sys.argv[1:]))
print 'Ryckaert-Bell coefficients:\n',C

A=array([
        1,   0, 1./2,    0, 3./8,     0,
        0,  -1,    0,-3./4,    0, -5./8,
        0,   0, 1./2,    0, 1./2,     0,
        0,   0,    0,-1./4,    0,-5./16,
        0,   0,    0,    0, 1./8,     0,
        0,   0,    0,    0,    0,-1./16,
        ]).reshape((6,6))

F=dot(A,C)
print 'Fourier coefficients:\n',F

x = arange(0,2*pi,0.01)
VF = F[0] + F[1]*cos(x) + F[2]*cos(2*x) + F[3]*cos(3*x) + \
    F[4]*cos(4*x) + F[5]*cos(5*x)
VRB = C[0] - C[1]*cos(x) + C[2]*(cos(x))**2 - C[3]*(cos(x)**3) + \
    C[4]*(cos(x))**4 - C[5]*(cos(x))**5

plot(x,VF,'b-')
plot(x,VRB,'g-')

ylabel('energy (kcal/mol)')
xlabel('torsional angle (rad)')

savefig('foo.pdf')
