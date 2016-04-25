% Check that we have a constrained state 
function checkConstrainedStateT(name,value)
    % Set the error message
    err = sprintf( ...
        'The %s argument must have type Constrained.State.t.', ...
        name);

    % Check the values
    try
        checkEqualityConstrainedStateT(name,value);
        checkInequalityConstrainedStateT(name,value);
    catch
        error(err);
    end
end
