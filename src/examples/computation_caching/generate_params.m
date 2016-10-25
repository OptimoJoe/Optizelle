% Generates all the parameters necessary for the simple convection diffusion
% optimization example
function params = generate_params()
    % Set the random seed
    randn('seed',2);

    % Set the problem size
    params.nx = 1000;

    % Set our exact solution 
    f = 4;
    m = 10;
    params.u = @(x)cos(2*pi*f*x)+m*x; 
    params.x_star = [ 1.23 ; 4.56];

    % Find our forcing function, x1 u'' + x2 u'
    params.f_analytic = @(x) ...
        params.x_star(1)*(-(2*pi*f)^2*cos(2*pi*f*x)) ...
        + params.x_star(2)*(-(2*pi*f)*sin(2*pi*f*x)+m);

    % Find our discretized domain
    params.a = 0;
    params.b = 1;
    params.dx = (params.b-params.a)/(params.nx+1);
    params.omega = linspace(params.a+params.dx,params.b-params.dx,params.nx)';

    % Discretieze our rhs
    params.f = params.f_analytic(params.omega);

    % Set our dirichlet boundary conditions
    dirichlet = params.u([params.a;params.b]);

    % Grab the adustment to the rhs from the boundary
    params.Ahat = generate_diffusion_bd(params.nx,params.dx,dirichlet);
    params.Bhat = generate_convection_bd(params.nx,params.dx,dirichlet);

    % Grab the differential operators
    params.A = generate_diffusion(params.nx,params.dx);
    params.B = generate_convection(params.nx,params.dx);

    % Grab some data based on our known solution and add some noise
    params.d = params.u(params.omega);
    pnoise = 0.25;
    noise = randn(params.nx,1);
    noise = (norm(params.d)/norm(noise))*noise;
    params.d = params.d + pnoise*noise;

    % Generate some indexing functions for the full-space problem
    params.idx.k=(1:2)';
    params.idx.u=(3:params.nx+2)';

    % Whether or not we're using the approximate Schur preconditioner for
    % the full-space problem
    params.approx_schur = 0;
end
