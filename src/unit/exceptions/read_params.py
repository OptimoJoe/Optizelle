# Tests error handling when reading invalid parameters 

import Optizelle 
import numpy

# Create some type shortcuts
XX = Optizelle.Rm

# Set the parameter file name 
fname = "bad_params.json"

# Allocate memory for an initial guess
x = numpy.array([1.2,2.3])

# Create an optimization state
state = Optizelle.Unconstrained.State.t(Optizelle.Rm,x)

msg = ""
#---Exception0---
# Read parameters from file
try:
    Optizelle.json.Unconstrained.read(XX,fname,state);
except Optizelle.Exception.t as e:
    # Convert the error message to a string 
    msg = e.message

    # Print the error message directly
    print e
#---Exception1---

# If we don't throw an exception above, throw an error
if len(msg)==0:
    raise Optizelle.Exception.t("Error catching missed our bad parameter")
