% Creates the 1-D finite difference operator for convection
function A = generate_convection(nx,dx)
    A = spdiags([-ones(nx,1) ones(nx,1)]/dx,-1:0,nx,nx);
end
