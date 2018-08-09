% Checks that an input is in a valid enumerated range
function checkEnumRange(name,enum,value)
    % Grab the field names in the enumerated type
    names = fieldnames(enum);

    % Grab the values out of the enumerated type using the field names.  We
    % grab one less because it's probably a function like toString.
    values = {};
    for i=1:length(names)-1
        values{i} = getfield(enum,names{i});
    end

    % Now, check if the number is in this list
    values = cell2mat(values);
    if sum(ismember(values,value))==0
        error(sprintf('The %s member is outside the valid enumated range.',...
            name));
    end
end
