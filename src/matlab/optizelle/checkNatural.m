% Checks that an input is a natural number 
function checkNatural(name,value)
    if ~isinteger(value) || value < 0
        error(sprintf('The %s member must be a natural number.',name));
    end
end
