#!/usr/bin/env python

from MagicUniverse import *
import sys

# preprocess_parallel_mesh(4); sys.exit(1)

NU = 0.1
RHO = 1.0
NSTEP = 60000
NDIAG = 100 # 1000
NVTKFREQ = 10
NCONFFREQ = -1
NDUMPFREQ = 10000
U_BAR = 0.01
Re = 5.0 #10.0
Pe = 10.0
C0 = 0.01
C1 = 1.0
UNFREEZE_TIME = 50000
GROWTIME = 2000 #1000
RESTART = False

R = Re * NU / 2 / U_BAR
DIFFUSIVITY = U_BAR * R / Pe

MagicBegins()

# define universe

u = Universe()
s = Scale()
m = Mesh()
f = Fluid()
c = Fluid()
a = Atom()
t = Tracker()

u.addItems([s,m,f,c,a,t])

u.setTitle('Re=5_RBC=5')
u.setNumberOfSteps(NSTEP)
u.setTemperature(0.0)
u.setStateRestart(RESTART) 
u.setStateDumpFrequency(NDUMPFREQ)

u.create()

# set params

s.set(name='MonoScale', mesh=m,  actors=[f,c,a,t])

m.setRegularMesh(True)

m.setPeriodicity('001')
m.setDomainDecomposition(7)
m.setPartitionAlongXYZ(2,2,256)


t.setDiagnosticFrequency(NDIAG)
# t.setDataShow(velocity=True)
# t.setMapDirections('zx')

t.setVtkDump(True, meshtype='unstructured', frequency=NVTKFREQ, start=50000, stop=54000)

f.setName('BloodFlow')
f.setCollisionType('BGK')
f.setDensityUniform(RHO)
f.setViscosity(NU)
f.setHomogeneousForce( 0., 0., 1.e-4 )
f.setFluidOnParticle(True)
f.setParticleOnFluid(True)

# NEWLY ADDED
c.setName('Bolus')
c.setDensityUniform(1.0)
c.setCollisionType('BGK')
c.setViscosity(DIFFUSIVITY) #changed from setDiffusivity
c.setAdvector(f)
c.setADR(True)
c.setFluidOnParticle(False)
c.setParticleOnFluid(True)

a.setName('RBC')
a.setIdentifier('molecule')
a.setConfigurationFile('atom.inp')
a.setLoop(outercycle=1, innercycle=1, innertimestep=1.)
a.setEnvironmentCouplingMethod('delta particle')
a.setConfigurationDump(start=0,frequency=NCONFFREQ)
a.setRotationalMotion(True)
# a.setForceAllPairs(True)

f.setInletOutletMethod('closed')
c.setInletOutletMethod('closed')
# f.setInletOutletMethod('equilibrium')
# f.setInletOutletMethod('closed')

f.setStabilizeLB(True)

u.decorate()

#f.setIOValue('inlet', 1, U_BAR)
#f.setIOValue('outlet', 2, U_BAR)
#c.setIOValue('inlet', 1, C0)
#c.setIOValue('outlet', 2, C0)

nx, ny, nz = m.getBox()
nx = int(nx); ny = int(ny); nz = int(nz)
profile2 = c.getArray(nx*ny*nz)
myid = get_myproc()

if not RESTART:
    profile2 = c.getArray(nx*ny*nz)
    for k in range(1,nz+1):
        for j in range(1,ny+1):
            for i in range(1,nx+1):

                ifl = m.getLocator(i,j,k)

                if k > nz/8. - 3 and k < nz/8. + 3: # create BOLUS
                    profile2[ifl] = C1
                else:
                    profile2[ifl] = C0

    c.setDensityProfile(profile2)

    # freeze the fluid and drug until released
    f.setFreeze(True)
    c.setFreeze(True)

    # cap RBC forces to a large roof to avoid instabilities
    # a.setCapForces(True, forcecap=1.e4, torquecap=1.e4, velcap=0.4, angvelcap=0.6)
    a.setCapForces(True, forcecap=1.e2, torquecap=1.e2, velcap=0.4, angvelcap=0.1)
    # set a robust friction coefficient for initial equilibration of RBC
    a.setGamma(gammaT=0.1, gammaR=0.1)
    a.setZeroVelocity()

for itime in u.cycle():

  if not RESTART:

    # gently increase excluded volume
    if itime <= GROWTIME:
        tscale = 0.5 * (1. + float(itime)/float(GROWTIME)) 
        if myid==0 and itime%10==0: 
            print 'Rescaling interaction to', tscale
            sys.stdout.flush()
        a.scaleVdwParameters(tscale)

    #  increase plasma-RBC coupling
    elif itime == int(1.5*GROWTIME):
        if myid==0:
            print 'Unfreezing fluid'
        f.setFreeze(False)
        a.setGamma(gammaT=0.0001, gammaR=0.0001)
        a.setZeroVelocity()

    # RBC free to move with right coupling
    elif itime == int(2.*GROWTIME):
        if myid==0:
            print 'RBC at right coupling'
        a.setGamma(gammaT=0.001, gammaR=0.001)
        a.setZeroVelocity()

    # RBC free to move with right coupling
    elif itime == int(2.5*GROWTIME):
        if myid==0:
            print 'RBC at right coupling'
        a.setGamma(gammaT=0.005, gammaR=0.005)
        a.setZeroVelocity()

    elif itime == int(UNFREEZE_TIME):
        if myid==0:
            print 'unfreezing drug'
        c.setFreeze(False)
  else:
    # allow the drug to move 
    if itime > int(2.5*GROWTIME):
        if myid==0:
            print 'RBC at right coupling'
        a.setGamma(gammaT=0.005, gammaR=0.005)
        a.setZeroVelocity()

    if itime == int(UNFREEZE_TIME):
        if myid==0:
            print 'unfreezing drug'
        c.setFreeze(False)

  # BC: stripe of given drug density
  c.setDensityStripe('z', nz, C0) 

  u.animate()

MagicEnds()