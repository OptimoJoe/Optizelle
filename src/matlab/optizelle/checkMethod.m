% Check that we have a method in a structure 
function result = checkMethod(name,value)
    result = isfield(value,name) && isa(getfield(value,name),'function_handle');
end
