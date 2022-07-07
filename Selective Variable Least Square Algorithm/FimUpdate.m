function [W_new b_new] = FimUpdate(W_old,b_old,a_new,input_new)
% Output: Updated FIM and Control Input Vector
% FIM Update
W_new = W_old + a_new'*a_new;
b_new = b_old + a_new'*input_new;
end

