% This tests our ability to capture and release from the optimization state
function unconstrained()

%---Import0---
global Optizelle
setupOptizelle();
%---Import1---

% Create some type shortcuts
XX = Optizelle.Rm;
msg = Optizelle.Messaging;

serialize=@(x)'';
deserialize=@(x,y)[];
check=@(x)0;
%---Serialization0---
Optizelle.json.Serialization.serialize('register',serialize,check); 
Optizelle.json.Serialization.deserialize('register',deserialize,check);
%---Serialization1---
    
% Create some arbitrary vector in R^2
x = [1.2;2.3];
x0 = [2.3;1.2];

% Create a state based on this vector
%---State0---
state = Optizelle.Unconstrained.State.t(XX,msg,x);
%---State1---

% Read in some parameters
fname = 'blank.json';
%---ReadJson0--- 
state = Optizelle.json.Unconstrained.read(XX,msg,fname,state);
%---ReadJson1--- 
   
% Create a bundle of functions
%---Functions0---
fns = Optizelle.Unconstrained.Functions.t;
%---Functions1---

% Do a null optimization
%---Solver0---
state = Optizelle.Unconstrained.Algorithms.getMin(XX,msg,fns,state);
%---Solver1---

% Do a null optimization with a state manipulator 
smanip = Optizelle.StateManipulator;
%---SmanipSolver0---
state = Optizelle.Unconstrained.Algorithms.getMin( ...
    XX,msg,fns,state,smanip);
%---SmanipSolver1---
    
% Read and write the state to file
fname = 'restart.json';
%---WriteReadRestart0---
Optizelle.json.Unconstrained.write_restart( ...
    XX,msg,fname,state);
state = Optizelle.json.Unconstrained.read_restart( ...
    XX,msg,fname,x);
%---WriteReadRestart1---

% Do a release 
%---Release0---
xs = Optizelle.Unconstrained.Restart.X_Vectors;
reals = Optizelle.Unconstrained.Restart.Reals;
nats = Optizelle.Unconstrained.Restart.Naturals;
params = Optizelle.Unconstrained.Restart.Params;
[xs reals nats params] = Optizelle.Unconstrained.Restart.release( ...
    XX,msg,state);
%---Release1---

% Check that we have the correct number of vectors
if length(xs) ~= 6
    msg.error('The list xs contains the wrong number of vectors.');
end

% Modify some vectors 
xs{1}{2}=x0;

% Capture the state
%---Capture0---
state = Optizelle.Unconstrained.Restart.capture( ...
    XX,msg,state,xs,reals,nats,params);
%---Capture1---

% Check the relative error between the vector created above and the one
% left in the state
residual = XX.init(x);
residual = XX.copy(x0);
residual = XX.axpy(-1.,state.x,residual);
err=(sqrt(XX.innr(residual,residual)) ...
    /(1+sqrt(XX.innr(x0,x0))));

if err >= 1e-15
    msg.error('Too much error in the captured x');
end
