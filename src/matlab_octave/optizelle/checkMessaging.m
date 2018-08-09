% Check that we have a messaging object
function checkMessaging(name,value)
    if ~isa(value,'function_handle')
        error(sprintf('The %s member must be a Messaging object.',name));
    end
end
