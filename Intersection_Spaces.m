function Intersection_Basis = Intersection_Spaces(L_1,L_2,tol)


% this function computes the Intersection between two linear spaces L_1 and
% L_2 given by column basis vectors

Intersection_Basis = null([null(L_1',tol)'; null(L_2',tol)'],tol);
Size_B = size(Intersection_Basis);
for row = 1:Size_B(1)
    for column = 1:Size_B(2)
        if abs(Intersection_Basis(row,column)) < tol
            Intersection_Basis(row,column) = 0;
        end
    end
end



end