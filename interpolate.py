import math,copy
from numpy import *

def frange(xmin,xmax,npoints):
    if npoints > 0:
        d = float(xmax - xmin)/(npoints - 1)
    else:
        d = 0
    r = []
    for i in range(npoints):
        r.append(xmin + i*d)
    return r

def sgn(x):
    if x < 0:
        return -1
    else:
        return 1

class PiecewisePolynomial:
    def __init__(self):
        """xc is a sorted list of lists [x0 c0 c1 c2 ..]
        where x0 is the start of each segment and ci is
        the coefficient of the term (x - x0)^i.
        The whole function is assumed to be continous."""
        self.xc = []

    def start(self):
        return self.xc[0][0]

    def end(self):
        return self.xc[-1][0]

    def __call__(self,x):
        """Return the interpolated value at x. If x is outside of the given
        data set None is returned."""
        if x < self.start() or x > self.end():
            return None
        else:
            # do a binary search for the right interval
            i = 0
            j = len(self.xc)
            while 1:
                k = int((i + j)/2)
                if x < self.xc[k][0]:
                    j = k
                if x >= self.xc[k][0]:
                    i = k
                if j <= i+1:
                    break
            #evaluate the local interpolant
            dx = x - self.xc[i][0]
            val = 0
            for j in range(len(self.xc[i])-1,0,-1):
                val = val*dx + self.xc[i][j]
            return val
        
    def derivative(self):
        """Return a new piecewise polynomial that is the derivative of
        this one. Works only if the derivative is representable as a
        polynomial of order one less than self."""
        p = PiecewisePolynomial()
        for x in self.xc:
            xp = [x[0]] + x[2:]
            for i in range(2,len(xp)):
                xp[i] *= i
            p.xc.append(xp)
        return p

    def integral(self):
        """Return the integral of this function."""
        raise NotImplementedError()

    def zeros(self, epsilon = 1e-15):
        """Return a list of all isolated zeros of the function, except
        at the boundaries. Only implemented for orders 1-3."""
        roots = []
        for i in range(len(self.xc)-1):
            poly = self.xc[i][1:]
            L = self.xc[i+1][0] - self.xc[i][0]
            if len(poly) == 2:
                # linear polynomial
                if -poly[0]*sgn(poly[1]) > 0 and -poly[0]*sgn(poly[1]) <= L*abs(poly[1]):
                    roots.append(-poly[0]/poly[1] + self.xc[i][0])
            elif len(poly) == 3 or (len(poly) == 4 and abs(poly[3]) < epsilon):
                # quadratic
                p = poly[1]**2 - 4*poly[0]*poly[2]
                if p >= 0:
                    p = math.sqrt(p)
                    if sgn(poly[2])*(-poly[1] + p) > 0 and sgn(poly[2])*(-poly[1] + p) <= L*2*abs(poly[2]):
                        roots.append((-poly[1] + p)/(2*poly[2]) + self.xc[i][0])
                    if sgn(poly[2])*(-poly[1] - p) > 0 and sgn(poly[2])*(-poly[1] - p) <= L*2*abs(poly[2]):
                        roots.append((-poly[1] - p)/(2*poly[2]) + self.xc[i][0])
            elif len(poly) == 4:
                # cubic, assumes poly[3] != 0.
                C = poly[0] / poly[3]
                B = poly[1] / poly[3]
                A = poly[2] / poly[3]
                Q = (3*B - A**2)/9
                R = (9*A*B - 27*C -2*A**3)/54
                D = Q**3 + R**2
                if D > 0:
                    #one real root
                    K = R - math.sqrt(D)
                    if abs(K) < epsilon:
                        x = -A/3 + (R + math.sqrt(D))**(1/3.0)
                    elif K > 0:
                        x = -A/3 + (R + math.sqrt(D))**(1/3.0) + K**(1/3.0)
                    else:
                        x = -1
                    if x > 0 and x <= L:
                        roots.append(x + self.xc[i][0])
                else:
                    # three real roots
                    theta = math.acos(R/math.sqrt(-Q**3))
                    M = 2*math.sqrt(-Q)
                    x = M*math.cos(theta/3) - A/3
                    if x > 0 and x <= L:
                        roots.append(x + self.xc[i][0])
                    x = M*math.cos((theta + 2*math.pi)/3) - A/3
                    if x > 0 and x <= L:
                        roots.append(x + self.xc[i][0])
                    x = M*math.cos((theta + 4*math.pi)/3) - A/3
                    if x > 0 and x <= L:
                        roots.append(x + self.xc[i][0])
            else:
                raise NotImplementedError()
        roots.sort()
        return roots

    def addConstant(self, y):
        for p in self.xc:
            p[1] += y

def Cubic(xy):
    p = PiecewisePolynomial()
    """Create a cublic spline interpolation object using the points
        in the xy list"""
    xy = map(lambda p: [float(p[0]),float(p[1])],xy)
    xy.sort()
    n = len(xy)
    b = [0]*n
    c = [0]*n
    d = [0]*n
    if n > 2:
        d[0] = xy[1][0] - xy[0][0]
        c[1] = (xy[1][1] - xy[0][1])/d[0]
        for i in range(1,n-1):
            d[i] = xy[i+1][0] - xy[i][0]
            b[i] = 2*(d[i-1] + d[i])
            c[i+1] = (xy[i+1][1] - xy[i][1])/d[i]
            c[i] = c[i+1] - c[i]
        b[0] = -d[0]
        b[-1] = -d[-2]
        c[0] = 0
        c[-1] = 0
        if n > 3:
            c[0] = c[2]/(xy[3][0] - xy[1][0]) - c[1]/(xy[2][0] - xy[0][0])
            c[-1] = c[-2]/(xy[-1][0] - xy[-3][0]) - c[-3]/(xy[-2][0] - xy[-4][0])
            c[0] = c[0]*d[0]**2/(xy[3][0] - xy[0][0])
            c[-1] = -c[-1]*d[-2]**2/(xy[-1][0] - xy[-4][0])
        for i in range(1,n):
            t = d[i-1]/b[i-1]
            b[i] = b[i] - t*d[i-1]
            c[i] = c[i] - t*c[i-1]
        c[-1] = c[-1]/b[-1]
        for ib in range(2,n+1):
            i = n - ib
            c[i] = (c[i] - d[i]*c[i+1])/b[i]
        b[-1] = (xy[-1][1] - xy[-2][1])/d[-2] + d[-2]*(c[-2] + 2*c[-1])
        for i in range(0,n-1):
            b[i] = (xy[i+1][1] - xy[i][1])/d[i] - d[i]*(c[i+1] + 2*c[i])
            d[i] = (c[i+1] - c[i])/d[i]
            c[i] = 3*c[i]
        c[-1] = 3*c[-1]
        d[-1] = d[-2]
    else:
        b[0] = (xy[1][1] - xy[0][1])/(xy[1][0] - xy[0][0])
        c[0] = 0
        d[0] = 0
        b[1] = b[0]
        c[1] = 0
        d[1] = 0
    for i in range(len(xy)):
        p.xc.append(xy[i] + [b[i], c[i],  d[i]])
    return p

def spline(x,y,npoints):
    xy = array([x.tolist(), y.tolist()]).transpose()
    f = Cubic(xy)
    xi=[]; yi=[]
    for t in frange(f.start(),f.end(),npoints):
        xi.append(t)
        yi.append(f(t))
    return array(xi), array(yi)
