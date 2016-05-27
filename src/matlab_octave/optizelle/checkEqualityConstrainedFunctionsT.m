% Check that we have an equality constrained function bundle 
function checkEqualityConstrainedFunctionsT(name,value)
    % Set the error message
    err = sprintf( ...
        'The %s argument must have type EqualityConstrained.Functions.t.', ...
        name);

    % Check the unconstrained values
    try
        checkUnconstrainedFunctionsT(name,value);
    catch
        error(err);
    end

    % Check for the equality constrained values 
    if ~(checkFields({ ...
        'g', ...
        'PSchur_left', ...
        'PSchur_right'}, ...
        value))
        error(err);
    end
end
