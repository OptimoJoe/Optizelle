% A function that has free reign to manipulate or analyze the state.  This 
% should be used cautiously.
function smanip = getStateManipulator()
    smanip=struct('eval',@(fns,state,loc)state);
