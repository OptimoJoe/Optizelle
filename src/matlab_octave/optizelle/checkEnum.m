% Checks that an input is an enumerated type
function checkEnum(name,value)
    if ~isinteger(value) || value < 0
        error(sprintf('The %s member must be an enumerated type (natural.)',...
            name));
    end
end
