n = 20;

fun = @objfungrad;

x0 = 0.4*ones(n^2,1);
% x0 = [zeros(0.3*n,n); ones(0.4*n,n); zeros(0.3*n,n)]; x0 = x0(:);
% x0 = x;
Aeq = ones(1,n^2);
beq = .4*n^2;
lb = zeros(n^2,1);
ub = ones(n^2,1);
options = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',true,...
    'MaxIterations',1e6,...
    'StepTolerance',1e-16,...
    'ConstraintTolerance',1e-16);
x = fmincon(fun,x0,[],[],Aeq,beq,lb,ub,[],options);

figure; heatmap( reshape(x,n,n));

function [f,gradf] = objfungrad(v)
p = 3;
[T, K, Q] = grid_discretisatie_khoek(v,p);
f = cost(T, K);
% Gradient of the objective function:
if nargout  > 1 || 1 == 1
    gradf = afkookt(v,p);
end
end

function c = cost(T, K)
    pow = 16;
    c = sum(T(:).^pow)^(1/pow);
%     c = Q'*T/2.532752327646418e+03;
%     c = .5*T'*K'*T;
end