% Check that we have several fields in a structure
function result = checkFields(names,value)
    result = 1;
    for i=1:length(names)
        result = result && isfield(value,names{i});
    end
end
