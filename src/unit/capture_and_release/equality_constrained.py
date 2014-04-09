# This tests our ability to capture and release from the optimization state

import Optizelle 
#---Import0---
import Optizelle.EqualityConstrained.State
import Optizelle.EqualityConstrained.Functions
import Optizelle.EqualityConstrained.Algorithms
import Optizelle.EqualityConstrained.Restart
import Optizelle.json.EqualityConstrained
#---Import1---
import numpy
import math

# Create some type shortcuts
XX = Optizelle.Rm
YY = Optizelle.Rm
msg = Optizelle.Messaging()
    
# Create some arbitrary vector in R^2
x = numpy.array([1.2,2.3])
x0 = numpy.array([2.3,1.2])

# Create a different arbitrary vector in R^3
y = numpy.array([3.4,4.5,5.6])
y0 = numpy.array([5.6,4.5,3.4])

# Create a state based on this vector
state=Optizelle.EqualityConstrained.State.t(XX,YY,msg,x,y)

# Read and write the state to file
fname = "restart.json"
#---WriteReadRestart0---
Optizelle.json.EqualityConstrained.write_restart(
    XX,YY,msg,fname,state);
Optizelle.json.EqualityConstrained.read_restart(
    XX,YY,msg,fname,x,y,state);
#---WriteReadRestart1---

# Do a release 
#---Release0---
xs = Optizelle.EqualityConstrained.Restart.X_Vectors()
ys = Optizelle.EqualityConstrained.Restart.Y_Vectors()
reals = Optizelle.EqualityConstrained.Restart.Reals()
nats = Optizelle.EqualityConstrained.Restart.Naturals()
params = Optizelle.EqualityConstrained.Restart.Params()
Optizelle.EqualityConstrained.Restart.release(
    XX,YY,msg,state,xs,ys,reals,nats,params)
#---Release1---

# Check that we have the correct number of vectors
if len(xs) != 14:
    msg.error("The list xs contains the wrong number of vectors.")
if len(ys) != 5:
    msg.error("The list ys contains the wrong number of vectors.")

# Modify some vectors 
xs[0]=(xs[0][0],x0)
ys[0]=(ys[0][0],y0)

# Capture the state
#---Capture0---
Optizelle.EqualityConstrained.Restart.capture(
    XX,YY,msg,state,xs,ys,reals,nats,params)
#---Capture1---

# Check the relative error between the vector created above and the one
# left in the state
residual = XX.init(x)
XX.copy(x0,residual)
XX.axpy(-1.,state.x,residual)
err=(math.sqrt(XX.innr(residual,residual))
    /(1+math.sqrt(XX.innr(x0,x0))))
if err >= 1e-15:
    msg.error("Too much error in the captured x")

residual = YY.init(y)
YY.copy(y0,residual)
YY.axpy(-1.,state.y,residual)
err=(math.sqrt(YY.innr(residual,residual))
    /(1+math.sqrt(YY.innr(y0,y0))))
if err >= 1e-15:
    msg.error("Too much error in the captured y")
