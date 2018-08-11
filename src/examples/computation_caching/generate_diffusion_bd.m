% Creates the 1-D finite difference boundary term for diffusion.  We use this
% to adjust the RHS.
function Ahat = generate_diffusion_bd(nx,dx,dirichlet)
    Ahat = zeros(nx,1);
    Ahat(1) = dirichlet(1)/dx^2;
    Ahat(nx) = dirichlet(2)/dx^2;
end
