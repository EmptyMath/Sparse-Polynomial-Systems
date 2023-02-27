function B = Adding_A(Monomials,A,Degree,Variables)

% this function adds A to all Monomials

Number_A = size(A);
Number_Monomials = size(Monomials);

for i = 1:Number_Monomials(1)
    for k = 1:Number_A
        new_monomial = Monomials(i,1:Variables) + A(k,:);
      
        if ismember(new_monomial, Monomials(1:Number_Monomials(1),1:Variables),'rows') == 0

           
            Monomials = [Monomials ; new_monomial Degree];
            Number_Monomials = size(Monomials);
        end
    end
end

B = Monomials;

end