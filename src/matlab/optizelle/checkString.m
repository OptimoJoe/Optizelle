% Checks that we have a string object 
function checkString(name,value)
    if ~ischar(value)
        error(sprintf('The %s member must be a string.',name));
    end
end
