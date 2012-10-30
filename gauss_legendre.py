#!/usr/bin/env python

"""
gauss_legendre: perform N-point Gauss-Legendre quadrature integration
                to evaluate C6 dispersion coefficient. Input to
                routines are 2 lists with the isotropic polarizability
                alpha(iw) evaluated at N (N=8,10,12) frequency
                points. The discrete frequencies are specific to the
                respective quadrature.

available routines:
                gauss_legendre_8pt(A,B)
                gauss_legendre_10pt(A,B)
                gauss_legendre_12pt(A,B)
"""

def gauss_legendre_8pt(A,B):
    "Gauss-Legendre 8 point quadrature on interval [-1,1]"
    if len(A)!=8:
        return "Incorrect numer of data points for A vector"
    if len(B)!=8:
        return "Incorrect numer of data points for B vector"
    from math import pi
    w0=0.3
    weights=[
        0.10122853629,
        0.222381034453,  
        0.313706645878,
        0.362683783378,
        0.362683783378,
        0.313706645878,
        0.222381034453,
        0.10122853629
        ]
    abscissas=[
        -0.960289856498, 
        -0.796666477414, 
        -0.525532409916,
        -0.183434642496,  
         0.183434642496,  
         0.525532409916,
         0.796666477414,  
         0.960289856498 
        ]
    c6=0.0
    for i in range(len(A)):
        x=-abscissas[i]
        w=weights[i]
        a=A[i]
        b=B[i]
        c6 = c6 + w*a*b/(1 + x)**2
    c6=c6*6.0*w0/pi
    print "    c6 =", c6
    return

def gauss_legendre_10pt(A,B):
    "Gauss-Legendre 10 point quadrature on interval [-1,1]"
    if len(A)!=10:
        return "Incorrect numer of data points for A vector"
    if len(B)!=10:
        return "Incorrect numer of data points for B vector"
    from math import pi
    w0=0.3
    weights=[
        0.066671344308,   
        0.149451349151,  
        0.219086362516,
        0.26926671931,    
        0.295524224715,  
        0.295524224715,
        0.26926671931,    
        0.219086362516,  
        0.149451349151,
        0.066671344308
        ]
    abscissas=[
        -0.973906528517, 
         -0.865063366689, 
         -0.679409568299,
         -0.433395394129, 
         -0.148874338982,  
         0.148874338982,
         0.433395394129,  
         0.679409568299,  
         0.865063366689,
         0.973906528517 
         ]
    c6=0.0
    for i in range(len(A)):
        x=-abscissas[i]
        w=weights[i]
        a=A[i]
        b=B[i]
        c6 = c6 + w*a*b/(1 + x)**2
    c6=c6*6.0*w0/pi
    print "    c6 =", c6
    return

def gauss_legendre_12pt(A,B):
    "Gauss-Legendre 12 point quadrature on interval [-1,1]"
    if len(A)!=12:
        return "Incorrect numer of data points for A vector"
    if len(B)!=12:
        return "Incorrect numer of data points for B vector"
    from math import pi
    w0=0.3
    weights=[
        0.0471753363866,
        0.106939325995,
        0.160078328543,
        0.203167426723,
        0.233492536538,
        0.249147045813,
        0.249147045813,
        0.233492536538,
        0.203167426723,
        0.160078328543,
        0.106939325995,
        0.0471753363866
        ]
    abscissas=[
        -0.981560634247,
        -0.90411725637,
        -0.769902674194,
        -0.587317954287,
        -0.367831498998,
        -0.125233408511,
        0.125233408511,
        0.367831498998,
        0.587317954287,
        0.769902674194,
        0.90411725637,
        0.981560634247]
    c6=0.0
    for i in range(len(A)):
        x=-abscissas[i]
        w=weights[i]
        a=A[i]
        b=B[i]
        c6 = c6 + w*a*b/(1 + x)**2
    c6=c6*6.0*w0/pi
    print "    c6 =", c6
    return

if __name__ == '__main__':
    print __doc__
