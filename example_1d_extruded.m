% Interpolate on a simple 1-D problem 

% Grab the Optizelle library
global Optizelle;
setupOptizelle ();

% Set the name of the parameter file
fname = 'example_1d_extruded.json';

% Set the size of the problem 
ninput = 2;
nhidden = 4;
nsamples = 20;

% Generate some random data along a direction
dir = rand(ninput,1);
dir = dir/norm(dir);
x = rand(ninput,nsamples);
true_fn = @(x)cos(5*x'*dir);
y = zeros(1,nsamples);
for j = 1:nsamples
    y(j) = true_fn(x(:,j)); 
end

% Allocate memory for an initial guess
xx = randn(nhidden+nhidden*ninput+nhidden,1);

% Create an optimization state
state=Optizelle.Unconstrained.State.t( ...
    Optizelle.Rm,Optizelle.Messaging,xx);

% Read the parameters from file
state=Optizelle.json.Unconstrained.read(Optizelle.Rm,Optizelle.Messaging,...
    fname,state);

% Grab some lenses
lens=generate_lenses(ninput,nhidden);

% Generate the objective function
fns=Optizelle.Unconstrained.Functions.t;
fns.f = generate_objective(generate_hyperbolic(),lens,x,y);

% Interpolate our data 
state=Optizelle.Unconstrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Messaging,fns,state);

% Generate an interpolatory function based on this
ff = generate_interpolant(generate_hyperbolic(),lens,state.x);

% Plot the result
dir_plot = rand(ninput,1);
dir_plot = dir_plot/norm(dir_plot,2);
x_uniform=0:.01:1;
figure(1);
plot(x_uniform,true_fn(dir_plot*x_uniform),'*', ...
     x_uniform,ff.eval(dir_plot*x_uniform),'x');
legend('True','Interpolation')

% Plot a surface plot in 2-D
if ninput==2
    % Figure out the points we're evaluating on
    [X Y]=meshgrid(linspace(0,1,100),linspace(0,1,100));

    % Evaluate our functions at these points
    Z_true = zeros(100,100);
    Z_interp = zeros(100,100);
    x_surf=zeros(2,10000);
    for i = 1 : 100
        for j = 1 : 100
            Z_true(i,j)=true_fn([X(i,j);Y(i,j)]);
            Z_interp(i,j)=ff.eval([X(i,j);Y(i,j)]);
        end
    end

    % Plot the two results
    figure(2);
    surf(X,Y,Z_true);
    title('True');
    axis([0 1 0 1 -2 2]);
    
    figure(3);
    surf(X,Y,Z_interp);
    title('Interpolation');
    axis([0 1 0 1 -2 2]);
end

% Calculate the direction the MLP chose
A = lens.A.get(state.x);
dir_mlp = zeros(size(A));
for i=1:nhidden
    dir_mlp(i,:)=A(i,:)/norm(A(i,:));
end

% Plot the directions if in 2-D
if ninput==2
    figure(4);
    compass(dir_mlp(:,1),dir_mlp(:,2),'b');
    hold on;
    compass(dir(1),dir(2),'r')
    hold off;
    title('Directions of Information Determined by the MLP')
end

