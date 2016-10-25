% Scales a vector so that its min and max values lie between the specified
% values
function f = generate_scaling_operator(from_min,from_max,to_min,to_max)
    f.eval = @(x)scaling_eval(from_min,from_max,to_min,to_max,x);
    f.p = @(x,dx)scaling_p(from_min,from_max,to_min,to_max,dx);
    f.ps = @(x,dy)scaling_p(from_min,from_max,to_min,to_max,dy);
    f.pps = @(x,dx,dy)zeros(size(x));
end

% Evaluation
function result = scaling_eval(from_min,from_max,to_min,to_max,x)
    % Figure out the amount of data and dimensions
    [ndim ndata] = size(x);

    % Loop over the dimensions and normalize
    result=zeros(size(x));
    for i=1:ndim
       result(i,:) = ...
           (to_max(i)-to_min(i)) / (from_max(i)-from_min(i)) ...
           * (x(i,:)-from_min(i)) ...
           + to_min(i);
    end
end

% Derivative 
function result = scaling_p(from_min,from_max,to_min,to_max,dx)
    % Figure out the amount of data and dimensions
    [ndim ndata] = size(dx);

    % Loop over the dimensions and normalize
    result=zeros(size(dx));
    for i=1:ndim
       result(i,:) =(to_max(i)-to_min(i)) / (from_max(i)-from_min(i))*(dx(i,:));
    end
end
