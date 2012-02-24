% Returns four functions that correspond to the residual function,
% its derivative, the adjoint of its derivative, and the adjoint of
% its second derivative after being given the first direction.

function r=gen_residual(d)

    % Basic residual
    f=@(i,y)y-d{i};

    % Derivative of the residual
    fp=@(i,y,dy)dy;

    % Adjoint of the derivative of the residual
    fps=@(i,y,dz)dz;

    % Adjoint of the second derivative of the residual in the direction dy1
    % applied in the direction dy2
    fpps=@(i,y,dy1,dy2)sparse([],[],[],size(dy1,1),1);

    % Build the entire residual operator
    r.f=f;
    r.fp=fp;
    r.fps=fps;
    r.fpps=fpps;
    r.max_index=length(d);
