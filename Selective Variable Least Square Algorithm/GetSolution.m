function solution = GetSolution(V_obs,D_obs,A,b)
%GETSOLUTION 이 함수의 요약 설명 위치
%   자세한 설명 위치
c = inv(D_obs)*V_obs'*A'*b;
solution = V_obs*c;

end

