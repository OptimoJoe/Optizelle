% Check that we have an inequality constrained function bundle 
function checkInequalityConstrainedFunctionsT(name,value)
    % Set the error message
    err = sprintf( ...
        'The %s argument must have type InequalityConstrained.Functions.t.', ...
        name);

    % Check the unconstrained values
    try
        checkUnconstrainedFunctionsT(name,value);
    catch
        error(err);
    end

    % Check for the inequality constrained values 
    if ~(checkFields({ ...
        'h'}, ...
        value))
        error(err);
    end
end
