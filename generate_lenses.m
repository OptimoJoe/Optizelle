% Generates the lenses based on the size of the network
function lens=generate_lenses(ninput,nhidden)

    % Determine the size of each variable
    idx = struct( ...
        'alpha',nhidden,
        'A',nhidden*ninput,
        'b',nhidden);
    
    % Generate some indexing values
    idx = generate_raw_indexing(idx);

    % Weights on the output layer  
    lens.alpha.get=@(x)x(idx.alpha);
    lens.alpha.set=@(dx,x)update(x,dx,idx.alpha);

    % Weights on the hidden layer
    lens.A.get=@(x)reshape(x(idx.A),nhidden,ninput);
    lens.A.set=@(dx,x)update(x,dx(:),idx.A);

    % Bias on the hidden layer
    lens.b.get=@(x)x(idx.b);
    lens.b.set=@(dx,x)update(x,dx,idx.b);
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
