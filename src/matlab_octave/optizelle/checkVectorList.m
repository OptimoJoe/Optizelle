% Checks that an input is a cell array
function checkVectorList(name,value)
    if ~iscell(value)
        error(sprintf('The %s member must be a cell array of vectors.',name));
    end
end
