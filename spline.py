#!/usr/bin/env python
import interpolate, sys

def frange(xmin,xmax,npoints):
    if npoints > 0:
        d = float(xmax - xmin)/(npoints - 1)
    else:
        d = 0
    r = []
    for i in range(npoints):
        r.append(xmin + i*d)
    return r


def eval_and_print(f, xmin, xmax, npoints):
    for x in frange(xmin,xmax,npoints):
        print x, f(x)


def readxy(f):
    xy = []
    for l in f.readlines():
        try:
            n = map(float,l.split())
            if len(n) == 2:
                xy.append(n)
        except Exception:
            pass
    return xy

def main():
    if len(sys.argv) == 1:
        print "Usage: spline [OPTIONS] OPERATION [FILE]"
        print "Reads xy data from stdin, do cubic spline interpolation"
        print "and perform some operation. Write the result to stdout."
        print ""
        print "OPERATION is one of"
        print " eval    Evaluate the spline function (see OPTIONS)"
        print " evald   Evaluate the derivative of the spline function"
        print " evaldd  Evaluate the second derivative of the spline function"
        print " zeros   Calculate all isolated zeros of the spline function inside"
        print "         its range, i.e. not at the start and end points."
        print " extrema Calculate all points with f' = 0, print x, f(x), f''(x)"
        print " opt YTHRES Return a new set of points so that a future spline"
        print "         interpolation of those points will have a maximum absolute"
        print "         vertical error YTHRES, with respect to the initial set."
        print ""
        print "OPTIONS can be"
        print " -o FILE     Write output to FILE (default: standard output)"
        print " -xmin XMIN  Start evaluation at XMIN (used for spline eval)"
        print " -xmax XMAX  End evaluation at XMAX (used for spline eval)"
        print " -points N   Evaluate the function at N points (used for spline eval)"
        print ""
        print "By <ulfek@ifm.liu.se>"
    elif sys.argv[1] == "eval":
        xy = readxy(sys.stdin)
        f = interpolate.Cubic(xy)        
        eval_and_print(f,f.start(),f.end(),1000)
    elif sys.argv[1] == "evald":
        xy = readxy(sys.stdin)
        f = interpolate.Cubic(xy)
        fp = f.derivative()
        eval_and_print(fp,fp.start(),fp.end(),1000)
    elif sys.argv[1] == "evaldd":
        xy = readxy(sys.stdin)
        f = interpolate.Cubic(xy)
        fp = f.derivative().derivative()
        eval_and_print(fp,fp.start(),fp.end(),1000)
    elif sys.argv[1] == "zeros":
        xy = readxy(sys.stdin)
        f = interpolate.Cubic(xy)
        for x in f.zeros():
            print x
    elif sys.argv[1] == "extrema":
        xy = readxy(sys.stdin)
        f = interpolate.Cubic(xy)
        fp = f.derivative()
        fb = fp.derivative()
        for x in fp.zeros():
            print x,f(x),fb(x)
    else:
        print "Unknown options",sys.argv

main()
    
