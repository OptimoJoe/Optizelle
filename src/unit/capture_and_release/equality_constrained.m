% This tests our ability to capture and release from the optimization state
function equality_constrained()

global Optizelle
setupOptizelle();

% Create some type shortcuts
XX = Optizelle.Rm;
YY = Optizelle.Rm;
msg = Optizelle.Messaging;
    
% Create some arbitrary vector in R^2
x = [1.2;2.3];
x0 = [2.3;1.2];

% Create a different arbitrary vector in R^3
y = [3.4;4.5;5.6];
y0 = [5.6;4.5;3.4];

% Create a state based on this vector
state=Optizelle.EqualityConstrained.State.t(XX,YY,msg,x,y);
    
% Read and write the state to file
fname = 'restart.json';
%---WriteReadRestart0---
Optizelle.json.EqualityConstrained.write_restart( ...
    XX,YY,msg,fname,state);
state = Optizelle.json.EqualityConstrained.read_restart( ...
    XX,YY,msg,fname,x,y);
%---WriteReadRestart1---

% Do a release 
%---Release0---
xs = Optizelle.EqualityConstrained.Restart.X_Vectors;
ys = Optizelle.EqualityConstrained.Restart.Y_Vectors;
reals = Optizelle.EqualityConstrained.Restart.Reals;
nats = Optizelle.EqualityConstrained.Restart.Naturals;
params = Optizelle.EqualityConstrained.Restart.Params;
[xs ys reals nats params] = Optizelle.EqualityConstrained.Restart.release( ...
    XX,YY,msg,state,xs,ys,reals,nats,params);
%---Release1---

% Check that we have the correct number of vectors
if length(xs) ~= 14
    msg.error('The list xs contains the wrong number of vectors.');
end
if length(ys) ~= 5
    msg.error('The list ys contains the wrong number of vectors.');
end

% Modify some vectors 
xs{1}{2}=x0;
ys{1}{2}=y0;

% Capture the state
%---Capture0---
state = Optizelle.EqualityConstrained.Restart.capture( ...
    XX,YY,msg,state,xs,ys,reals,nats,params);
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
    msg.error('Too much error in the captured y')
end
