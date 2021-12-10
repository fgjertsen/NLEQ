% Test script to solve equations by Broyden update vs. Newtons method

% Eq:
%	x1^2	+	x2^2	+	x3^2	= 3
%	x1^2	+	x2^2	-	x3		= 1
%	x1		+	x2		+	x3		= 3

% Jac.
%	2x	2y	2z
%	2x	2y	-1
%	1	1	1

function nleq_example_comparison
	clc; close all;
	f = @(x, par) fx(x);
	J = @(x, par) jac(x);	

	% Initial guess for solution
	%x0 = [1 0 1]';
	x0 = [0.5 1.5 0.5]';
	n = length(x0);
	
	% Settings
	tol = 1.e-12;
	nItMax = 150;
	par = [];
	
	% Solve
	tic; fprintf('Doing Newton iterations...\t');
		xm_newton  = nleq_newton( x0, f, par, J, tol, nItMax);	toc;
		fprintf('Number of iterations needed:\t');
		nIter_newton  = size(xm_newton,  2) - 1;
		fprintf('%d\n', nIter_newton);
	tic; fprintf('Doing Broyden iterations...\t');
		xm_broyden = nleq_broyden(x0, f, par, J, tol, nItMax);	toc;
		fprintf('Number of iterations needed:\t');
		nIter_broyden  = size(xm_broyden,  2) - 1;
		fprintf('%d\n', nIter_broyden);
	
	% Plots
	h = figure;
	set( h, 'Position', [100 100 1000 500] );
	for i = 1:n
		subplot(1,n,i); plot(0:nIter_newton,	xm_newton(i,:), 'Color', 'blue'); hold on;
						plot(0:nIter_broyden,	xm_broyden(i,:), 'Color', 'red');
						plot(nIter_newton, 	xm_newton(i, end), 'Color', 'blue', 'Marker', '*');
						plot(nIter_broyden, xm_broyden(i, end), 'Color', 'red', 'Marker', '*');
			title(['x' num2str(i)]);
			legend('Newton', 'Broyden', 'Location', 'northeast');
	end

end

% Jacobian
function ret = jac(x)
ret = [ 2*x(1), 2*x(2), 2*x(3)	;...
		2*x(1),	2*x(2), -1		;...
		1,		1,		1		 ...
		];
end

% Equation
function ret = fx(x)
ret = [	x(1)^2 + x(2)^2 + x(3)^2 - 3	; ...
		x(1)^2 + x(2)^2 - x(3)   - 1	; ...
		x(1)   + x(2)   + x(3)   - 3	  ...
		];	
end