function Basis = Construct_Basis(B,New_Basis_Monomials,C_0,R,Variables,Monomials,tol)

% this function constructs a Basis B connected to 1 such that B + C_0 = R
% and B cap C_0 = 0, we have that B ans C_0 are subspaces of R

% we start with a basis B and only add new basis elements by multiplying by
% x_i, hence we get a basis connected to 1

B = [B];

New_Basis_Monomials_2 = [];


if Same_Spaces([B C_0],R,tol) == 0 
    Prolongation_Monomials = InverseCoefficient(New_Basis_Monomials,Monomials);
    Number_of_Monomials = length(Prolongation_Monomials);
    for k = 1:Number_of_Monomials

        % we are going to multiply each new basis monomial by x_i

        for i = 1:Variables
            Monomial = Prolongation_Monomials{k};
            Monomial{1}{2}(i) = Prolongation_Monomials{k}{1}{2}(i) + 1;
            New_Basis_Monomial = Coefficients(Monomial,Monomials);
            

            % we need to check if the new basis monomial is indeed linear
            % independent and that it is containted in R, if all of that is
            % the case then we add it to the basis and also remember that
            % this monomial is a 'new' one



            if rank([B C_0],tol) < rank([B New_Basis_Monomial C_0],tol) & ismember(New_Basis_Monomial',R',"rows") == 1
               B = [B New_Basis_Monomial];
               New_Basis_Monomials_2 = [New_Basis_Monomials_2 New_Basis_Monomial];
            end
        end
    end
    Basis = Construct_Basis(B,New_Basis_Monomials_2,C_0,R,Variables,Monomials,tol);
else
    Basis = B;

end

end