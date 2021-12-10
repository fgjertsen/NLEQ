% FUNCTION NAME:
%   nleq_broyden
%
% DESCRIPTION:
%   Find the solution for a nonlinear (set of) equation(s)
%		on the form Ax = b, by iteration, using Broydens method
%
% INPUT:
%   xv0		- (double, vector) 				Initial guess for solution
%   f   	- (fhandle(double, vector))		Function to solve, f(x, par) = 0
%	par		- (double, struct)				System specific parameters
%	J 		- (fhandle(double, matrix)) 	Function to calculate Jacobian, df/dx.	(optional)
%	tol		- (double, scalar)				Tolerance for iteration convergence 	(optional)
%	nItMax	- (integer, scalar)				Max number of iterations 				(optional)
%
% OUTPUT:
%   xm	- (double, vector) 	Solution to equation
%			Dimension is [n, nIterations]
%
function xm = nleq_broyden(xv0, f, par, J, tol, nItMax)

	% Initiate function values
	xv = xv0;
	n  = length(xv);
	fr = feval(f, xv, par);
	
	% Initiate Broyden approximation to Jacobian
	if not( isempty ( J ) )
		Br = inv( feval( J, xv, par) );
	else
		Br = eye(n);
	end
	M = zeros(n);
	
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
	
	% Tuning parameter of Broydens method
	tau = 1;
	
	% Start iteration loop
	bConverged = 0;
	nIt = 0;
	while ~bConverged
		
		% Save iterates
		xm(:, nIt+1) = xv;
		
		% Check for convergence
		if norm(fr) < tol || nIt >= nItMax
			bConverged = 1;
			break;
		end
		
		% Do "Newton Iteration" step
		pr = -Br * fr;
		xv = xv + tau*pr;
		
		% Update Jacobian approximation, using Broydens formula
		dfr = feval(f, xv, par) - fr;
		fr  = fr + dfr;
		oyp = Br*dfr - pr;
		pB	= pr' * Br;
		for i=1:n
			for j=1:n
				M(i,j) = oyp(i) * pB(j);
			end
		end
		Br = Br - M./(pr'*Br*dfr);
		
		% Increment loop counter
		nIt = nIt + 1;
		
	end
	
	% Return the solutions from the iterations
	xm = xm(:, 1:nIt+1);

end