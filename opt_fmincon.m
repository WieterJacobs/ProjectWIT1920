n = 30;

fun = @objfungrad;

x0 = 0.4*ones(n^2,1);
Aeq = ones(1,n^2);
beq = .4*n^2;
lb = zeros(n^2,1);
ub = ones(n^2,1);
options = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',true,...
    'MaxIterations',1e6,...
    'StepTolerance',1e-16);
%     'ConstraintTolerance',1e-6);
x = fmincon(fun,x0,[],[],Aeq,beq,lb,ub,[],options);

heatmap( reshape(x,n,n));

function [f,gradf] = objfungrad(v)
p = 3;
f = max(grid_discretisatie_khoek(v,p));
% Gradient of the objective function:
if nargout  > 1 || 1 == 1
    gradf = afkookt(v,p);
end
end