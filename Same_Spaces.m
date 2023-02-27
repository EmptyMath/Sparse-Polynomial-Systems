function T = Same_Spaces(L_1,L_2,tol)

% this function returns 1 if and only if L_1 and L_2 span the same vector
% space, we use a rank condition for this.

% L_1 and L_2 are given by a basis consisting of column vectors

C = [L_1, L_2];

if rank(C,tol) == rank(L_1,tol) && rank(C,tol) == rank(L_2,tol)
    T = 1;
else
    T = 0;
end

end