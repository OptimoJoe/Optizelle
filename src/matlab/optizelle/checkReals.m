% Check that we have a list of restart reals 
function checkReals(name,value)
    if ~iscell(value)
        error(sprintf('The %s argument must be a cell array.',name));
    end
    for i=1:length(value)
        checkString(sprintf('%s{%d}{1}',name,i),value{i}{1});
        checkFloat(sprintf('%s{%d}{2}',name,i),value{i}{2});
    end
end
