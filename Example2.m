f_1 = {{1,[2 11]},{-1,[1 0]}};
f_2 = {{1,[3 0]},{-1,[0 4]}};
f_3 = {{1,[1 2]},{1,[0 2]},{-2,[0 0]}};
Polynomials = {f_1,f_2,f_3};

A = [0 0; 1 0; 0 1; 0 2; 0 4; 2 11; 3 0];
tol = 1e-3;
Variables = 2;

Roots = Moment_Method_A_standard(Polynomials,A,Variables,tol)