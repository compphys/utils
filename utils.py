#!/usr/bin/env python
import numpy
from math import exp

def lorentzian(x,y,gamma,initial,final,n=1000):
    step=(x[-1]-x[0])/n
    xi=numpy.arange(initial,final,step)
    yi=numpy.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(y)):
            yi[i] = yi[i] + y[k] * gamma**2 / ( (xi[i]-x[k])**2 + gamma**2 )
    return xi,yi

def gaussian(x,y,sigma,initial,final,n=1000):
    step=(x[-1]-x[0])/n
    xi=numpy.arange(initial,final,step)
    yi=numpy.zeros(len(xi))
    for i in range(len(xi)):
        for k in range(len(y)):
            yi[i] = yi[i] + y[k] * exp( -(xi[i]-x[k])**2 / sigma**2 )
    return xi,yi

def array_reverse_order(x):
    xtmp = numpy.zeros(len(x))
    for i in range(len(x)):
        xtmp[i]=x[-1-i]
    return xtmp

def idip(d):
    if d=='X' or d=='x':
        id=1
    elif d=='Y' or d=='y':
        id=2
    elif d=='Z' or d=='z':
        id=3
    else:
        id=-1
    return id

def find_string_in_file(f,s):
    p=-1
    while p==-1:
        r=f.readline()
        p=r.find(s)
    return r

def skiplines(f,n):
    for i in range(n):
        f.readline()

def from_file(file):
    f = open(file, 'r')
    s = f.read()
    f.close()
    return s

def to_file(s,file):
    f = open(file, 'w')
    f.write(s)
    f.close()

def lines_from_file(file):
    return from_file(file).split('\n')

def load(file):
    from numpy import array
    f = open(file, 'r')
    a = []
    for l in f.readlines():
        try:
            n = map(float,l.replace(',',' ').split())
            if len(n)>0:
                a.append(n)
        except ValueError:
            pass
    return array(a)

def index_limits(x,lower,upper):
    i=0
    while x[i] < lower:
        i=i+1
    ilower=i
    while x[i] < upper:
        i=i+1
    iupper=i
    return ilower,iupper
