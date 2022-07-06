function obs_index = CheckObsParam(V_obs)
%CHECKOBSPARAM 이 함수의 요약 설명 위치
%   자세한 설명 위치
base = eye(3);
obs_index = [];
[R C] = size(V_obs);

    for i = 1:C
        for j = 1:3
            a(j) = abs(base(:,j)'*V_obs(:,i));
        end
        obs_index = [obs_index find(a == max(a))];
    end    
end