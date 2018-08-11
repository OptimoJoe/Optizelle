% Check that we have several methods in a structure
function result = checkMethods(names,value)
    result = 1;
    for i=1:length(names)
        result = result && checkMethod(names{i},value);
    end
end
