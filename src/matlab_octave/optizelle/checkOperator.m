% Check that we have a linear operator
function checkOperator(name,value)

    % Check for the appropriate fields and make sure they're functions
    if ~(checkMethod('eval',value))
        error(sprintf('The %s member must be an Operator.',name));
    end
end
