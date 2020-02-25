n = 20;
xnew = 0.4*ones(n); xnew = xnew(:);
gnew = afkookt(xnew);
lambda = 100;

iter = 2000;
for i = 1:iter
    xold = xnew;
    gold = gnew;
    xnew = xold - lambda*gold;
    xnew = scale(xnew);
    gnew = afkookt(xnew);
    lambda = abs((xnew-xold)'*(gnew - gold));
    lambda = lambda/norm(gnew - gold)^2;
    
    if mod(i,50) == 0
        figure(1);
        v = reshape(xnew,n,n);
        v = v(2:end-1,2:end-1);
        heatmap([v fliplr(v); flipud(v) rot90(v,2)]);
        drawnow()
        figure(2);
        g = reshape(gnew,n,n);
        heatmap(g);
        drawnow()
    end
end


function x = scale(x)
x = max(x,0);
x = min(x,1);

x = x/mean(x)*0.4;
end