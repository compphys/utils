#!/usr/bin/env python

# Physical constants
from math import pi

a0    = 5.291772108e-11
m     = 9.1093826e-31
e     = 1.60217653e-19
hbar  = 1.05457168e-34
Eh    = 4.35974417e-18
pi    = 3.1415926535
c     = 2.99792458e8
epsilon0 = 8.8541878e-12
alpha = e*e/(4*pi*epsilon0*hbar*c)

if __name__ == '__main__':
    for label in dir():
        if label.find('__') == -1:
            print '%10s' % (label + '='),
            exec 'print ' + label
