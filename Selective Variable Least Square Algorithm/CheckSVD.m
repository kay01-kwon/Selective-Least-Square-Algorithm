function [V_obs D_obs] = CheckSVD(V,D)
%CHECKSVD 이 함수의 요약 설명 위치
%   자세한 설명 위치
SV_max = 50;
Lambda = [D(1,1);D(2,2);D(3,3)];
Lambda_max = max(Lambda);
SV1 = sqrt(abs(Lambda_max/Lambda(1)));
SV2 = sqrt(abs(Lambda_max/Lambda(2)));

Lambda_obs = [];
V_obs = [];

if abs(SV1) < SV_max
    V_obs = [V_obs V(:,1)];
    Lambda_obs = [Lambda_obs Lambda(1)];
end

if abs(SV2) < SV_max
    V_obs = [V_obs V(:,2)];
    Lambda_obs = [Lambda_obs Lambda(2)];
end

V_obs = [V_obs V(:,3)];
Lambda_obs = [Lambda_obs Lambda(3)];
D_obs = diag(Lambda_obs);
end

