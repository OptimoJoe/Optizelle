% Check that we have a valid Euclidean-Jordan algebra
function checkEuclidean(name,value)

    % Check for the appropriate fields and make sure they're functions
    if ~(checkMethods({'prod','id','linv','barr','srch','symm'))
        error(sprintf('The %s member must be a Euclidean-Jordan algebra.',...
            name));
    end
end
