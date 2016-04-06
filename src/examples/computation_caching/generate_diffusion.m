% Creates the 1-D finite difference operator for diffusion
function A = generate_diffusion(nx,dx)
    A = spdiags([ones(nx,1) -2*ones(nx,1) ones(nx,1)]/dx^2,-1:1,nx,nx);
end
