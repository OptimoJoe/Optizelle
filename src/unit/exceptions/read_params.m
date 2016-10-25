% Tests error handling when reading invalid parameters 
function read_params()
    % Grab the Optizelle library
    global Optizelle;
    setupOptizelle();

    % Create some type shortcuts
    XX = Optizelle.Rm;

    % Set the parameter file name 
    fname = 'bad_params.json';

    % Allocate memory for an initial guess
    x = [1.2,2.3];

    % Create an optimization state
    state=Optizelle.Unconstrained.State.t(Optizelle.Rm,x);

    % Try to catch the rror 
    msg = '';
    %---Exception0---
    try
        state = Optizelle.json.Unconstrained.read(XX,fname,state);
    catch e
        % Convert the error message into a string
        msg = e.message;

        % Print the message directly
        disp(e.message);
    end
    %---Exception1---

    % If we don't throw an exception above, throw an error
    if length(msg)==0
        error('Error catching missed our bad parameter')
    end
end
