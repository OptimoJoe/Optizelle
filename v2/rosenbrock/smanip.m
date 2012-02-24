function manipulated=smanip(original)
    % Plot the Rosenbrock function and the current iterate
    if original.iter==1
    	clf
	[X Y]=meshgrid(-1.5:.1:1.5,-1.5:.1:1.5);
	contour(X,Y,myfunc2(X,Y),100)  
	hold on
    end
    plot(original.u(1),original.u(2),'xr');

    % Don't manipulate the state
    manipulated=original;
