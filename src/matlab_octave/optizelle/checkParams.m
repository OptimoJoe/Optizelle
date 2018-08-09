% Check that we have a list of restart parameters
function checkParams(name,value)
    if ~iscell(value)
        error(sprintf('The %s argument must be a cell array.',name));
    end
    for i=1:length(value)
        checkString(sprintf('%s{%d}{1}',name,i),value{i}{1});
        checkString(sprintf('%s{%d}{2}',name,i),value{i}{2});
    end
end
