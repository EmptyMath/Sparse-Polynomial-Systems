function [At,b,c,K] = SDP_Input_A(Polynomials,Monomials_A,A,Variables,order)

% this function computes the SDP input in SeDuMi format for the moment
% method, the constant term gets a sign '-' 
% the monomials are given as row vectors
% the polynomials are given in the form Polynomials {f_1,...,f_m}, and each
% polynomial is given as f_i = {{[c_1],[x^alpha_1]},...,{[c_k],[x^alpha_k]}}
% Variables is an integer

% the monomials_A have to have the same degree as the order, not more!

% switch b and c???

N = size(Monomials_A);
Number_of_monomials = N(1);
Number_of_polynomials = length(Polynomials);
Index_of_1 = 1;
while isequal(Monomials_A(Index_of_1,1:Variables),zeros(1,Variables)) == 0
    Index_of_1 = Index_of_1 + 1;
end

b = zeros(1,Number_of_monomials);
b(Index_of_1) = 1;

% which size does c have?

% Writing the At, which consists of the matrices A_j represented as column
% vectors

% assuming Index_of_1 = 1

At = [1 zeros(1,Number_of_monomials-1)];
Degree_A = Degree_of_Polynomials_A(Polynomials,Variables,A);

for beta = 1:Number_of_monomials
    for i = 1:Number_of_polynomials
        
        
        A_alpha = [];
        if Monomials_A(beta,Variables+1) <= order - Degree_A(i)
        
           N = length(Polynomials{i});
           for alpha = 1:Number_of_monomials
               c_i_alpha_beta = 0;
               % checking if alpha - beta is an exponent of h_i
               
               for k = 1:N
                   
                   
                   if isequal(Monomials_A(alpha,1:Variables) - Monomials_A(beta,1:Variables), Polynomials{i}{k}{2}) == 1
                       c_i_alpha_beta = Polynomials{i}{k}{1};
                   end
               end

               A_alpha = [A_alpha, c_i_alpha_beta];
           end
        end
        
        rows = size(A_alpha);
        if rank(At) < rank([At; A_alpha])
        At = [At ; A_alpha];
        end
        
    end
       
    
    
    
   
end

% get the size for C out of At or A_alpha
size_At = size(At);
C = [1 zeros(1,size_At(1)-1)] ;



Monomials_of_order_t_2 = 0;
for i = 1:Number_of_monomials
    if Monomials_A(i,Variables+1) <= order/2
        Monomials_of_order_t_2 = Monomials_of_order_t_2 + 1;
    end
end


rows_At = size(At);
At = [At zeros(rows_At(1),Monomials_of_order_t_2^2)];


b = [b  zeros(1,Monomials_of_order_t_2 ^2)]';
C = [C  zeros(1,Monomials_of_order_t_2 ^2)]';


k = 0;
% Momonet matrix should be psd using new variables, and setting them as
% entries in the moment matrix by the constraint = 0
for alpha = 1:Number_of_monomials
if Monomials_A(alpha,Variables+1) <= order
    k_alpha = 0;
    Monomials_A(alpha,:);
    New_constrains = [];
    y_alpha = zeros(1,Number_of_monomials);
    y_alpha(1,alpha) = 1;

    for i = 1:Number_of_monomials
        T = 0;
        for j = 1:Number_of_monomials
            
            
            if Monomials_A(i,Variables+1) <= order/2 & Monomials_A(j,Variables + 1) <= order/2 & T == 0
                
                if isequal(Monomials_A(i,1:Variables) + Monomials_A(j,1:Variables), Monomials_A(alpha,1:Variables)) 
                    
                   
                    
                    
                    z_alpha = zeros(1,Monomials_of_order_t_2^2);
                    z_alpha(1,Monomials_of_order_t_2*(i-1)+j) = -1;
                    New_constrains = [New_constrains ; y_alpha z_alpha];
                    T = 1;
                end
            end
        end
    end
%if rank(At) < rank([At ; New_constrains])
       At = [At ; New_constrains];
       size(At);
%end
   


end

end



L = [At C];
L = L(any(L,2),:);
size_L = size(L);
At = L(1:size_L(1),1:size_L(2)-1);
c = L(1:size_L(1),size_L(2));

% declare the free and sdp variables

K.f = Number_of_monomials;
if order == 1
    K.l = 1;
    K.s = 0;
else
    K.l = 0;
    K.s = Monomials_of_order_t_2;
end

end