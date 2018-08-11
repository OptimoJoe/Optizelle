% Check that we have a valid vector space
function checkVectorSpace(name,value)

    % Check for the appropriate fields and make sure they're functions
    if ~(checkMethods({'init','copy','scal','zero','axpy','innr','rand'},value))
        error(sprintf('The %s member must be a vector space.',name));
    end
end
