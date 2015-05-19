% Check that we have a scalar-valued function 
function checkScalarValuedFunction(name,value)

    % Check for the appropriate fields and make sure they're functions
    if ~(checkMethods({'eval','grad','hess_vec'},value)) 
        error(sprintf('The %s member must be a ScalarValuedFunction.',name));
    end
end
