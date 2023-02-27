function Roots = Eigenvalue_Method(B,M,Multiplication_Maps,Variables)

% this function computes the roots with the border basis approach.
% B, the border basis is given as row vectors and M is a matrix.

% Compute the left eigenvectors of the Multiplication matrix, therefore
% transpose the matrix. By transposing V, we get the left eigenvectors

[V,D] = eig(M');
Left_Eigenvectors = V';

% we are now looking where the monomial 1 stands in the basis B

Index_of_1 = 1;

while isequal(B(Index_of_1,:),zeros(1,Variables)) == 0
    Index_of_1 = Index_of_1 + 1;
end

Cardinality_of_B = size(B);

% computing the roots as explained in the thesis

Roots = [];

for k = 1:Cardinality_of_B(1)
    root = zeros(1,Variables);
    Eigenvector = Left_Eigenvectors(k,:);
    Constant = Eigenvector(Index_of_1);
    
    for i = 1:Variables
        x_i = zeros(1,Variables);
        x_i(i) = 1;        

        if ismember(x_i,B,'rows') == 1

            % we need the index of x_i
            Index_of_x = 1;
            while isequal(B(Index_of_x,:),x_i) == 0
                Index_of_x = Index_of_x + 1;
            end
            root(i) = Eigenvector(Index_of_x) / Constant;
        else

            % making use of eq (4) in the thesis

% Are we looking at the right normal form????
            
            Normal_Form_x_i = Multiplication_Maps(:,(i-1)*Cardinality_of_B(1) + 1);

            Eig = Eigenvector / Constant;
            for j = 1:Cardinality_of_B(1)
                root(i) = root(i) + Normal_Form_x_i(j) * Eig(j);
            end

        end
        
    end
    Roots = [Roots; root];

end


end