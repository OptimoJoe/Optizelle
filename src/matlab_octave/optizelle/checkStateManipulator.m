% Check that we have a state manipulator
function checkStateManipulator(name,value)

    % Check for the appropriate fields and make sure they're functions
    if ~(checkMethod('eval',value))
        error(sprintf('The %s member must be a StateManipulator object.',name));
    end
end
