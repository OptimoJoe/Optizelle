% Check that we have a constrained function bundle
function checkConstrainedFunctionsT(name,value)
    % Set the error message
    err = sprintf( ...
        'The %s argument must have type Constrained.Functions.t.', ...
        name);

    % Check the values
    try
        checkEqualityConstrainedFunctionsT(name,value);
        checkInequalityConstrainedFunctionsT(name,value);
    catch
        error(err);
    end
end
