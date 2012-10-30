#!/usr/bin/env python

# Conversion constants

au2ev = 27.2114
ev2cm = 8065.540901
au2cm = au2ev*ev2cm
amu2au = 1822.888479031408
au2debye = 2.54177
au2kcal = 627.5095
au2ang=0.5291772083
au2s=2.418884326505e-17

if __name__ == '__main__':
    for label in dir():
        if label.find('__') == -1:
            print '%10s' % (label + '='),
            exec 'print ' + label
