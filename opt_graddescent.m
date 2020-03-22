n = 30;
p = 3;
xnew = 0.4*ones(n); xnew = xnew(:);
gnew = afkookt(xnew,p);
lambda = 100;

iter = 10000;
for i = 1:iter
    xold = xnew;
    gold = gnew;
    xnew = xold - lambda*gold;
    xnew = scale(xnew);
    gnew = afkookt(xnew,p);
%     gnew = fin_dif(xnew,p);
    lambda = abs((xnew-xold)'*(gnew - gold));
    lambda = lambda/norm(gnew - gold)^2;
    
    if mod(i,50) == 0
        figure(1);
        v = reshape(xnew,n,n);
        v = v(2:end-1,2:end-1);
        heatmap(v);
        heatmap([v fliplr(v); flipud(v) rot90(v,2)]);
        drawnow()
%         figure(2);
%         g = reshape(gnew,n,n);
%         heatmap(g);
%         drawnow()
%         figure(3);
%         T = grid_discretisatie_khoek(xnew,p);
%         heatmap(reshape(T,n+1,n+1));
    end
end



function x = scale(x)
x = x/mean(x)*0.4;
x = max(x,0);
x = min(x,1);
x = x/mean(x)*0.4;
end

% function g = filter(g,n)
% g = reshape(g,n,n);
% for i = 2:n-1
%     for j = 2:n-1
%         M1 = g(i-1:i+1,j-1:j+1);
%         M2 = [1 2 1; 2 4 2; 1 2 1];
%         g(i,j) = g(i,j)+sum(sum(M1.*M2));
%         g(i,j) = g(i,j)/16;
%     end
% end
% g = g(:);
% end