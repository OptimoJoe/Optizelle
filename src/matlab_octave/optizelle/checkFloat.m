% Checks that an input is a floating-point number
function checkFloat(name,value)
    if ~isreal(value)
        error(sprintf('The %s member must be a floating point.',name));
    end
end
