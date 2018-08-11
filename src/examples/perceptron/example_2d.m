% Interpolate on a simple 2-D problem

% Grab the Optizelle library
global Optizelle;
setupOptizelle ();

% Set the name of the parameter file
fname = 'example_2d.json';

% Set whether we're plotting or not
do_plot = 0;

% Set the size of the problem
ninput = 2;
nhidden = 7;
nsamples = 50;

% Determine if we black out a quadrant (no data)
blackout = 0;
x_black = [0.5;1];
y_black = [0.5;1];

% Generate some random data for the product of some sines and cosines along
% skewed axes
dir = rand(ninput,2);
for j=1:2
    dir(:,j) = dir(:,j)/norm(dir(:,j));
end

% Figure out where we're sampling
x = rand(ninput,nsamples);
if blackout
    for j=1:nsamples
        while x(1,j) >= x_black(1) && x(1,j) <= x_black(2) && ...
              x(2,j) >= y_black(1) && x(2,j) <= y_black(2)
            x(:,j)=rand(ninput,1);
        end
    end
end

true_fn = @(x)cos(4*x'*dir(:,1)).*sin(2*x'*dir(:,2));
y = zeros(1,nsamples);
for j = 1:nsamples
    y(j) = true_fn(x(:,j));
end

% Generate scalings based on this data
scaling = generate_scaling(x,y);

% Grab some lenses
[lens idx]=generate_lenses(ninput,nhidden);

% Allocate memory for an initial guess
xx = randn(idx.size,1);

% Create an optimization state
state=Optizelle.Unconstrained.State.t(Optizelle.Rm,xx);

% Read the parameters from file
state=Optizelle.json.Unconstrained.read(Optizelle.Rm,fname,state);

% Generate the objective function
fns=Optizelle.Unconstrained.Functions.t;
fns.f = generate_objective(generate_hyperbolic(),lens,x,y,scaling);

% Interpolate our data
state=Optizelle.Unconstrained.Algorithms.getMin( ...
    Optizelle.Rm,Optizelle.Messaging.stdout,fns,state);

% Generate an interpolatory function based on this
ff = generate_interpolant(generate_hyperbolic(),lens,state.x,scaling);

% Return if we're not plotting
if ~do_plot
    return;
end

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

    % Create a small sphere
    [sx sy sz] = sphere(5);
    sx = 0.02*sx;
    sy = 0.02*sy;
    sz = 0.02*sz;

    % Plot the two results
    figure(2);
    h = surf(X,Y,Z_true);
    set(h, 'cdata',zeros(100))

    hold on
    h = surf(X,Y,Z_interp);
    set(h, 'cdata',0.5* ones(100))
    hold off
    title('True (blue) vs Interpolated (green)');
    axis([0 1 0 1 -2 2]);

    % Plot where the data is
    hold on
    for j = 1:nsamples
        h = surf(sx+x(1,j),sy+x(2,j),sz+y(j));
        set(h, 'cdata',ones(6))
    end
    hold off
end

% Calculate the direction the MLP chose
A = lens.A.get(state.x);
dir_mlp = zeros(size(A));
for i=1:nhidden
    dir_mlp(i,:)=A(i,:)/norm(A(i,:));
end

% Plot the directions if in 2-D
if ninput==2
    figure(3);
    compass(dir_mlp(:,1),dir_mlp(:,2),'b');
    hold on;
    compass(dir(:,1),dir(:,2),'r')
    hold off;
    title('Directions of Information Determined by the MLP')
end

