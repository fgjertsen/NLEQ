% FUNCTION NAME:
%   nleq_newton
%
% DESCRIPTION:
%   Find the solution for a nonlinear (set of) equation(s)
%		on the form Ax = b, by iteration, using Newtons method
%
% INPUT:
%   xv0		- (double, vector) 				Initial guess for solution
%   f   	- (fhandle(double, vector))		Function to solve, f(x, par) = 0
%	par		- (double, struct)				System specific parameters
%	J 		- (fhandle(double, matrix)) 	Function to calculate Jacobian, df/dx.
%	tol		- (double, scalar)				Tolerance for iteration convergence 	(optional)
%	nItMax	- (integer, scalar)				Max number of iterations 				(optional)
%
% OUTPUT:
%   xm	- (double, matrix) 	Solution to equation, at each iteration
%			Dimension is [n, nIterations]
%
function xm = nleq_newton(xv0, f, par, J, tol, nItMax)

	% Initiate function values
	xv = xv0;
	n  = length(xv);
	
	% Make sure tolerance is set
	if isempty(tol)
		tol = 1.e-10;
	end
	
	% Make sure max. iterations is set
	if isempty(nItMax)
		nItMax = 1000;
	end
	
	% Allocate solutions matrix, with sufficient space
	xm = zeros(n, nItMax+1);
	
	% Start iteration loop
	bConverged = 0;
	nIt = 0;
	while ~bConverged
		
		% Save iterates
		xm(:, nIt+1) = xv;
		
		% Evaluate nonlin. eq.		
		fr = feval( f, xv, par );
		
		% Check for convergence
		if norm(fr) < tol || nIt >= nItMax
			bConverged = 1;
			break;
		end
		
		% Do Newton Iteration step
		pr = -inv( feval( J, xv, par ) ) * fr;
		xv = xv + pr;
		
		% Increment loop counter
		nIt = nIt + 1;
		
	end
	
	% Return the solutions from the iterations
	xm = xm(:, 1:nIt+1);

end
