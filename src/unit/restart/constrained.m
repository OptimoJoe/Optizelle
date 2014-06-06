% This tests our ability to capture and release from the optimization state
function constrained()

global Optizelle
setupOptizelle();

% Create some type shortcuts
XX = Optizelle.Rm;
YY = Optizelle.Rm;
ZZ = Optizelle.Rm;
msg = Optizelle.Messaging;
    
% Create some arbitrary vector in R^2
x = [1.2;2.3];
x0 = [2.3;1.2];

% Create a different arbitrary vector in R^3
y = [3.4;4.5;5.6];
y0 = [5.6;4.5;3.4];

% Create a different arbitrary vector in R^4
z = [6.7;7.8;8.9;9.10];
z0 =[9.10;8.9;7.8;6.7];

% Create a state based on this vector
%---State0---
state = Optizelle.Constrained.State.t(XX,YY,ZZ,msg,x,y,z);
%---State1---

% Read in some parameters
fname = 'blank.json'
%---ReadJson0--- 
state = Optizelle.json.Constrained.read(XX,YY,ZZ,msg,fname,state);
%---ReadJson1--- 
   
% Create a bundle of functions
%---Functions0---
fns = Optizelle.Constrained.Functions.t;
%---Functions1---

% Do a null optimization
%---Solver0---
state = Optizelle.Constrained.Algorithms.getMin(XX,YY,ZZ,msg,fns,state);
%---Solver1---

% Do a null optimization with a state manipulator 
smanip = Optizelle.StateManipulator;
%---SmanipSolver0---
state = Optizelle.Constrained.Algorithms.getMin( ...
    XX,YY,ZZ,msg,fns,state,smanip);
%---SmanipSolver1---

% Read and write the state to file
fname = 'restart.json';
%---WriteReadRestart0---
Optizelle.json.Constrained.write_restart( ...
    XX,YY,ZZ,msg,fname,state);
state = Optizelle.json.Constrained.read_restart( ...
    XX,YY,ZZ,msg,fname,x,y,z);
%---WriteReadRestart1---

% Do a release 
%---Release0---
xs = Optizelle.Constrained.Restart.X_Vectors;
ys = Optizelle.Constrained.Restart.Y_Vectors;
zs = Optizelle.Constrained.Restart.Z_Vectors;
reals = Optizelle.Constrained.Restart.Reals;
nats = Optizelle.Constrained.Restart.Naturals;
params = Optizelle.Constrained.Restart.Params;
[xs ys zs reals nats params] = Optizelle.Constrained.Restart.release( ...
    XX,YY,ZZ,msg,state,xs,zs,reals,nats,params);
%---Release1---

% Check that we have the correct number of vectors
if length(xs) ~= 14 
    msg.error('The list xs contains the wrong number of vectors.');
end
if length(ys) ~= 5
    msg.error('The list ys contains the wrong number of vectors.');
end
if length(zs) ~= 3
    msg.error('The list zs contains the wrong number of vectors.');
end

% Modify some vectors 
xs{1}{2}=x0;
ys{1}{2}=y0;
zs{1}{2}=z0;

% Capture the state
%---Capture0---
state = Optizelle.Constrained.Restart.capture( ...
    XX,YY,ZZ,msg,state,xs,ys,zs,reals,nats,params);
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

residual = YY.init(y);
residual = YY.copy(y0);
residual = YY.axpy(-1.,state.y,residual);
err=(sqrt(YY.innr(residual,residual)) ...
    /(1+sqrt(YY.innr(y0,y0))));
if err >= 1e-15
    msg.error('Too much error in the captured y');
end

residual = ZZ.init(z);
residual = ZZ.copy(z0);
residual = ZZ.axpy(-1.,state.z,residual);
err=(sqrt(ZZ.innr(residual,residual)) ...
    /(1+sqrt(ZZ.innr(z0,z0))));
if err >= 1e-15
    msg.error('Too much error in the captured z')
end
