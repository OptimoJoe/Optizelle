# This tests our ability to capture and release from the optimization state

import Optizelle 
import numpy
import math

# Create some type shortcuts
XX = Optizelle.Rm
YY = Optizelle.Rm
ZZ = Optizelle.Rm
msg = Optizelle.Messaging.stdout
    
# Create some arbitrary vector in R^2
x = numpy.array([1.2,2.3])
x0 = numpy.array([2.3,1.2])

# Create a different arbitrary vector in R^3
y = numpy.array([3.4,4.5,5.6])
y0 = numpy.array([5.6,4.5,3.4])

# Create a different arbitrary vector in R^4
z = numpy.array([6.7,7.8,8.9,9.10])
z0 = numpy.array([9.10,8.9,7.8,6.7])

# Create a state based on this vector
#---State0---
state = Optizelle.Constrained.State.t(XX,YY,ZZ,x,y,z)
#---State1---

# Read in some parameters
fname = "blank.json"
#---ReadJson0--- 
Optizelle.json.Constrained.read(XX,YY,ZZ,fname,state)
#---ReadJson1--- 
   
# Create a bundle of functions
#---Functions0---
fns = Optizelle.Constrained.Functions.t()
#---Functions1---

# Do a null optimization
state.f_x = 1.0
#---Solver0---
Optizelle.Constrained.Algorithms.getMin(XX,YY,ZZ,msg,fns,state)
#---Solver1---

# Do a null optimization with a state manipulator 
smanip = Optizelle.StateManipulator()
#---SmanipSolver0---
Optizelle.Constrained.Algorithms.getMin(XX,YY,ZZ,msg,fns,state,smanip)
#---SmanipSolver1---

# Read and write the state to file
fname = "restart.json"
#---WriteReadRestart0---
Optizelle.json.Constrained.write_restart(XX,YY,ZZ,fname,state);
Optizelle.json.Constrained.read_restart(XX,YY,ZZ,fname,x,y,z,state);
#---WriteReadRestart1---

# Do a release 
#---Release0---
xs = Optizelle.Constrained.Restart.X_Vectors()
ys = Optizelle.Constrained.Restart.Y_Vectors()
zs = Optizelle.Constrained.Restart.Z_Vectors()
reals = Optizelle.Constrained.Restart.Reals()
nats = Optizelle.Constrained.Restart.Naturals()
params = Optizelle.Constrained.Restart.Params()
Optizelle.Constrained.Restart.release(
    XX,YY,ZZ,state,xs,ys,zs,reals,nats,params)
#---Release1---

# Check that we have the correct number of vectors
if len(xs) != 14:
    raise Optizelle.Exception.t(
        "List xs contains the wrong number of vectors")
if len(ys) != 5:
    raise Optizelle.Exception.t(
        "List ys contains the wrong number of vectors")
if len(zs) != 3:
    raise Optizelle.Exception.t(
        "List zs contains the wrong number of vectors")

# Modify some vectors 
xs[0]=(xs[0][0],x0)
ys[0]=(ys[0][0],y0)
zs[0]=(zs[0][0],z0)

# Capture the state
#---Capture0---
Optizelle.Constrained.Restart.capture(
    XX,YY,ZZ,state,xs,ys,zs,reals,nats,params)
#---Capture1---

# Check the relative error between the vector created above and the one
# left in the state
residual = XX.init(x)
XX.copy(x0,residual)
XX.axpy(-1.,state.x,residual)
err=(math.sqrt(XX.innr(residual,residual))
    /(1+math.sqrt(XX.innr(x0,x0))))
if err >= 1e-15:
    raise Optizelle.Exception.t("Too much error in the captured x")

residual = YY.init(y)
YY.copy(y0,residual)
YY.axpy(-1.,state.y,residual)
err=(math.sqrt(YY.innr(residual,residual))
    /(1+math.sqrt(YY.innr(y0,y0))))
if err >= 1e-15:
    msg("Too much error in the captured y")

residual = ZZ.init(z)
ZZ.copy(z0,residual)
ZZ.axpy(-1.,state.z,residual)
err=(math.sqrt(ZZ.innr(residual,residual))
    /(1+math.sqrt(ZZ.innr(z0,z0))))
if err >= 1e-15:
    msg("Too much error in the captured z")
