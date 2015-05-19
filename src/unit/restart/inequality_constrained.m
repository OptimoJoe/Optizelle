% This tests our ability to capture and release from the optimization state
function inequality_constrained()

global Optizelle
setupOptizelle();

% Create some type shortcuts
XX = Optizelle.Rm;
ZZ = Optizelle.Rm;
msg = Optizelle.Messaging;
    
% Create some arbitrary vector in R^2
x = [1.2;2.3];
x0 = [2.3;1.2];

% Create a different arbitrary vector in R^4
z = [6.7;7.8;8.9;9.10];
z0 =[9.10;8.9;7.8;6.7];

% Create a state based on this vector
%---State0---
state = Optizelle.InequalityConstrained.State.t(XX,ZZ,msg,x,z);
%---State1---

% Read in some parameters
fname = 'blank.json';
%---ReadJson0--- 
state = Optizelle.json.InequalityConstrained.read(XX,ZZ,msg,fname,state);
%---ReadJson1--- 
   
% Create a bundle of functions
%---Functions0---
fns = Optizelle.InequalityConstrained.Functions.t;
%---Functions1---

% Do a null optimization
state.f_x = 1.0;
%---Solver0---
state = Optizelle.InequalityConstrained.Algorithms.getMin(XX,ZZ,msg,fns,state);
%---Solver1---

% Do a null optimization with a state manipulator 
smanip = Optizelle.StateManipulator;
%---SmanipSolver0---
state = Optizelle.InequalityConstrained.Algorithms.getMin( ...
    XX,ZZ,msg,fns,state,smanip);
%---SmanipSolver1---
    
% Read and write the state to file
fname = 'restart.json';
%---WriteReadRestart0---
Optizelle.json.InequalityConstrained.write_restart( ...
    XX,ZZ,msg,fname,state);
state = Optizelle.json.InequalityConstrained.read_restart( ...
    XX,ZZ,msg,fname,x,z);
%---WriteReadRestart1---

% Do a release 
%---Release0---
xs = Optizelle.InequalityConstrained.Restart.X_Vectors;
zs = Optizelle.InequalityConstrained.Restart.Z_Vectors;
reals = Optizelle.InequalityConstrained.Restart.Reals;
nats = Optizelle.InequalityConstrained.Restart.Naturals;
params = Optizelle.InequalityConstrained.Restart.Params;
[xs zs reals nats params] = Optizelle.InequalityConstrained.Restart.release( ...
    XX,ZZ,msg,state);
%---Release1---

% Check that we have the correct number of vectors
if length(xs) ~= 6 
    msg.error('The list xs contains the wrong number of vectors.');
end
if length(zs) ~= 3
    msg.error('The list zs contains the wrong number of vectors.');
end

% Modify some vectors 
xs{1}{2}=x0;
zs{1}{2}=z0;

% Capture the state
%---Capture0---
state = Optizelle.InequalityConstrained.Restart.capture( ...
    XX,ZZ,msg,state,xs,zs,reals,nats,params);
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

residual = ZZ.init(z);
residual = ZZ.copy(z0);
residual = ZZ.axpy(-1.,state.z,residual);
err=(sqrt(ZZ.innr(residual,residual)) ...
    /(1+sqrt(ZZ.innr(z0,z0))));
if err >= 1e-15
    msg.error('Too much error in the captured z')
end
