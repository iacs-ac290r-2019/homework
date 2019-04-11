#!/usr/bin/env python

from muphy2wrapper import *
import sys

RHO = 0.22
NSTEP = 40000
NU = None
um = None
VIN = 0.4 # m/s : Hyperemia

def units():

    global NSTEP,RHO,NU,um

    um = UnitsMapping()
    um.setResolution(100) #vox/cm
    um.setPhysicalMassDensity(1000.) #Kg/m3
    um.setPhysicalViscosity(4.e-6) #m2/s
    um.setCharacteristicVelocity(0.2) #m/s
    um.setCharacteristicPressure_mmHg(100.) # mmHg

    um.setCharacteristicMach(0.025) # best practice for resolution 300
    # um.setCharacteristicMach(0.035)
    # um.setCharacteristicMach(0.045)

    # Hyperhemia Set Here
    um.addInlet(id=1, velocity=VIN, area=14.51) #mm2 and m/s

    um.addOutlet(id=2, name='OM1', area=2.39)
    um.addOutlet(id=3, name='D1',  area=1.77)
    um.addOutlet(id=4, name='LCX', area=2.08)
    um.addOutlet(id=5, name='LAD', area=2.23)
    um.addOutlet(id=6, name='OM2', area=2.33)
    um.addOutlet(id=7, name='X0',  area=1.48)
    um.addOutlet(id=8, name='X1',  area=2.30)
    um.addOutlet(id=9, name='X2',  area=1.07)

    um.setOutletMethodOutflow()

    um.workout()

    NU = um.getLatticeViscosity()

    print '\nSimulation Total Physical Time:',NSTEP*um.getTimeStep(),' s\n'

# check that *.ios header is ok: pressure for inlet, flow for outlets
def check(f):

    global NSTEP,RHO,NU,um

    for key in um.inlets.keys():
        id = um.inlets[key].Id
        bctype = f.getIOBCType('inlet',id)
        if bctype[:8] != 'pressure':
            print 'Inlet BC must be pressure!'
            sys.exit(1)

    for key in um.outlets.keys():
        id = um.outlets[key].Id
        bctype = f.getIOBCType('outlet',id)
        if bctype[:4] != 'flow':
            print 'Outlet BC must be flow!'
            sys.exit(1)
