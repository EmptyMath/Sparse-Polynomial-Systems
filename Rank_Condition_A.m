function [T,order] = Rank_Condition_A(D,t,lambda,monomials,tol,Variables,A)

% this function checks the rank conditions explained in Thm 2.46 and gives
% the integer s which fulfills one of these conditions

T = 0;
d = ceil(D/2);
order = 0;

% checking the rank condition ii) in Thm 2.46

if d <= t
    for s = d:t
        
        
        M_order_s = Moment_Matrix_A(lambda,s,monomials,Variables,A);
        
        
        M_order_s_d = Moment_Matrix_A(lambda,s-d,monomials,Variables,A);
        rank(M_order_s_d,tol);
        if rank(M_order_s_d,tol) == rank(M_order_s,tol)
            order = s;
            T = 1;
            return
        end
    end
end

% checking the rank condition i) in Thm 2.46

if D <= t
    for s = D:t
        
        
        M_order_s = Moment_Matrix_A(lambda,s,monomials,Variables,A);
        rank(M_order_s);
        
        M_order_s_1 = Moment_Matrix_A(lambda,s-1,monomials,Variables,A);
     
        if rank(M_order_s_1,tol) == rank(M_order_s,tol)
            order = s;
            
            T = 1;
            return
        end
    end
end

end