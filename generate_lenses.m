% Generates the lenses based on the size of the network
function [get set]=generate_lenses(nhidden,noutput)

    % Determine the size of each variable
    idx = struct( ...
        'alpha',noutput,
        'A',noutput*nhidden,
        'b',noutput);
    
    % Generate some indexing values
    idx = generate_raw_indexing(idx);

    % Create some functions that get the specified values 
    get.alpha=@(x)x(raw_idx.alpha);
    get.A=@(x)reshape(x(raw_idx.A),noutput,nhidden);
    get.b=@(x)x(raw_idx.b);

    % Create some functions that set the specified values 
    set.alpha=@(x,dx)update(x,dx,idx.alpha);
    set.A=@(x,dx)update(x,dx(:),idx.A);
    set.b=@(x,dx)update(x,dx,idx.b);
end

% Creates indices for elements specified in the given structure
function idx = generate_raw_indexing(idx)
    % Grab these field names
    fn = fieldnames(idx);

    % Loop over each of the field names while keeping track of the current 
    % starting index number
    curr = 1;
    for i=1:length(fn)
        % Get the size of the current field
        m = idx.(fn{i});
       
        % Overwrite the size of the field with its indexing
        idx.(fn{i}) = curr:(curr+m-1);

        % Update the current index 
        curr = curr + m;
    end

    % Add a special size element to the indexing function
    idx.size = curr-1;
end

% Does an update operation 
function x=update(x,dx,idx)
    x(idx)=dx;
end
