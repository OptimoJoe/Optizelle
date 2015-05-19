% Check that we have a messaging object 
function checkMessaging(name,value)

    % Check for the appropriate fields and make sure they're functions
    if ~(checkMethods({'print','error'},value)) 
        error(sprintf('The %s member must be a Messaging object.',name));
    end
end
