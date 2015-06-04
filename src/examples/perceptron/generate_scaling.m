% Generates scalings to bound the inputs between given values.  If the user
% doesn't specify bounds, we default to on the inputs [0,1] and [-1,1] on
% the outputs.
function scaling = generate_scaling(x,y,scaling)
    % Determine the scalings
    scaling.x.from.min = eval('scaling.x.from.min','min(x'')''');
    scaling.x.from.max = eval('scaling.x.from.max','max(x'')''');
    
    scaling.x.to.min = eval('scaling.x.to.min','zeros(size(x,1),1)');
    scaling.x.to.max = eval('scaling.x.to.max','ones(size(x,1),1)');
    
    scaling.y.from.min = eval('scaling.y.from.min','min(y'')''');
    scaling.y.from.max = eval('scaling.y.from.max','max(y'')''');
    
    scaling.y.to.min = eval('scaling.y.to.min','-ones(size(y,1),1)');
    scaling.y.to.max = eval('scaling.y.to.max','ones(size(y,1),1)');
end
