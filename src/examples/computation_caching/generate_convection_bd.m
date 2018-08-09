% Creates the 1-D finite difference boundary term for convection.  We use this
% to adjust the RHS.
function Ahat = generate_convection_bd(nx,dx,dirchlet)
    Ahat = zeros(nx,1);
    Ahat(1) = -dirchlet(1)/dx;
end
