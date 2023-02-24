f_1 = {{1,[0 0]},{-1,[-1 1]}};
f_2 = {{1,[2 0]},{-2,[0 0]}};
Polynomials = {f_1,f_2};
Variables = 2;

A = [0 0; 1 0; -1 0; 0 1; 0 -1];
tol = 1e-7;

Roots = Moment_Method_A(Polynomials,A,Variables,tol)