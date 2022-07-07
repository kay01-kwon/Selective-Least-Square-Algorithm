function solution = GetRecursiveSolution(V_obs,D_obs,b)
% Output: solution
c = inv(D_obs)*V_obs'*b;
solution = V_obs*c;
end

