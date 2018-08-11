% Check that we have a vector-valued function
function checkVectorValuedFunction(name,value)

    % Check for the appropriate fields and make sure they're functions
    if ~(checkMethods({'eval','p','ps','pps'},value))
        error(sprintf('The %s member must be a VectorValuedFunction.',name));
    end
end
