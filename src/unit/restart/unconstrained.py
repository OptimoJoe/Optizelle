# This tests our ability to capture and release from the optimization state

#---Import0---
import Optizelle 
#---Import1---

import numpy
import math

# Create some type shortcuts
XX = Optizelle.Rm
msg = Optizelle.Messaging.stdout

serialize=lambda x:""
deserialize=lambda x,y:[]
vector_type=type(serialize)
#---Serialization0---
Optizelle.json.Serialization.serialize.register(serialize,vector_type)
Optizelle.json.Serialization.deserialize.register(deserialize,vector_type)
#---Serialization1---
    
# Create some arbitrary vector in R^2
x = numpy.array([1.2,2.3])
x0 = numpy.array([2.3,1.2])

# Create an unconstrained state based on this vector
#---State0---
state = Optizelle.Unconstrained.State.t(XX,x)
#---State1---
    
# Read in some parameters
fname = "blank.json"
#---ReadJson0--- 
Optizelle.json.Unconstrained.read(XX,fname,state)
#---ReadJson1--- 
   
# Create a bundle of functions
#---Functions0---
fns = Optizelle.Unconstrained.Functions.t()
#---Functions1---

# Do a null optimization
state.f_x = 1.0
#---Solver0---
Optizelle.Unconstrained.Algorithms.getMin(XX,msg,fns,state)
#---Solver1---

# Do a null optimization with a state manipulator 
smanip = Optizelle.StateManipulator()
#---SmanipSolver0---
Optizelle.Unconstrained.Algorithms.getMin(XX,msg,fns,state,smanip)
#---SmanipSolver1---

# Read and write the state to file
fname = "restart.json"
#---WriteReadRestart0---
Optizelle.json.Unconstrained.write_restart(XX,fname,state);
Optizelle.json.Unconstrained.read_restart(XX,fname,x,state);
#---WriteReadRestart1---

# Do a release 
#---Release0---
xs = Optizelle.Unconstrained.Restart.X_Vectors()
reals = Optizelle.Unconstrained.Restart.Reals()
nats = Optizelle.Unconstrained.Restart.Naturals()
params = Optizelle.Unconstrained.Restart.Params()
Optizelle.Unconstrained.Restart.release(XX,state,xs,reals,nats,params)
#---Release1---

# Check that we have the correct number of vectors
if len(xs) != 6:
    raise Optizelle.Exception.t(
        "The list xs contains the wrong number of vectors.")

# Modify some vectors 
xs[0]=(xs[0][0],x0)

# Capture the state
#---Capture0---
Optizelle.Unconstrained.Restart.capture(XX,state,xs,reals,nats,params)
#---Capture1---

# Check the relative error between the vector created above and the one
# left in the state
residual = XX.init(x)
XX.copy(x0,residual)
XX.axpy(-1.,state.x,residual)
err=(math.sqrt(XX.innr(residual,residual))
    /(1+math.sqrt(XX.innr(x0,x0))))

if err >= 1e-15:
    msg("Too much error in the captured x")
