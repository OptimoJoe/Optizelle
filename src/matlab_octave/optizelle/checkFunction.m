% Checks that an input is a function
function checkFunction(name,value)
    if ~isa(value,'function_handle')
        error(sprintf('The %s member must be a function.',name));
    end
end
