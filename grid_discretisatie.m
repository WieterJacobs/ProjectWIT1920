function T = grid_discretisatie(v)
% F*T = Q
n = size(v,1);
deltax = 0.01;
deltan = deltax/n;
k = v*65 + (1 - v)*0.2;
F = zeros((n+2)^2);
Q = zeros(n+2); Q = Q(:);

for i = 1:n
    for j = 1:n
        Q(pos(i,j)) = 2*deltan^2*10^11;
        counter = 0;
        if i-1 ~= 0
            counter = counter + 1;
            F(pos(i,j),pos(i-1,j)) = -4*k(i-1,j);
            F(pos(i,j),pos(i-2,j)) = k(i-1,j);
            F(pos(i,j),pos(i,j)) = k(i-1,j);
            F(pos(i,j),pos(i-1,j-1)) = k(i-1,j);
            F(pos(i,j),pos(i-1,j+1)) = k(i-1,j);
        end
        if i+1 ~= n+1
            counter = counter + 1;
            F(pos(i,j),pos(i+1,j)) = -4*k(i+1,j);
            F(pos(i,j),pos(i,j)) = k(i+1,j);
            F(pos(i,j),pos(i+2,j)) = k(i+1,j);
            F(pos(i,j),pos(i+1,j-1)) = k(i+1,j);
            F(pos(i,j),pos(i+1,j+1)) = k(i+1,j);
        end
        if j-1 ~= 0
            counter = counter + 1;
            F(pos(i,j),pos(i,j-1)) = -4*k(i,j-1);
            F(pos(i,j),pos(i-1,j-1)) = k(i,j-1);
            F(pos(i,j),pos(i+1,j-1)) = k(i,j-1);
            F(pos(i,j),pos(i,j-2)) = k(i,j-1);
            F(pos(i,j),pos(i,j)) = k(i,j-1);
        end
        if j+1 ~= n+1
            counter = counter + 1;
            F(pos(i,j),pos(i,j+1)) = -4*k(i,j+1);
            F(pos(i,j),pos(i-1,j+1)) = k(i,j+1);
            F(pos(i,j),pos(i+1,j+1)) = k(i,j+1);
            F(pos(i,j),pos(i,j)) = k(i,j+1);
            F(pos(i,j),pos(i,j+2)) = k(i,j+1);
        end
        F(pos(i,j),pos(i,j)) = counter*4*k(i,j);
        F(pos(i,j),pos(i-1,j)) = -counter*k(i,j);
        F(pos(i,j),pos(i+1,j)) = -counter*k(i,j);
        F(pos(i,j),pos(i,j-1)) = -counter*k(i,j);
        F(pos(i,j),pos(i,j+1)) = -counter*k(i,j);
    end
end

Dir = 1+round(.3*n):n+2-round(.3*n);
for j = 0:n+1
    F(pos(0,j),pos(0,j)) = 1;
    F(pos(0,j),pos(1,j)) = -1;
%     Q(pos(0,j)) = 0;
    F(pos(n+1,j),pos(n+1,j)) = 1;
    F(pos(n+1,j),pos(n,j)) = -1;
%     Q(pos(n+1,j)) = 0;
end
for i = 1:n
    if ismember(i, Dir)
        F(pos(i,0),pos(i,0)) = 1;
        Q(pos(i,0)) = 293;
        F(pos(i,n+1),pos(i,n+1)) = 1;
        Q(pos(i,n+1)) = 293;
    else
        F(pos(i,0),pos(i,0)) = 1;
        F(pos(i,0),pos(i,1)) = -1;
%         Q(pos(i,0)) = 0;
        F(pos(i,n+1),pos(i,n+1)) = 1;
        F(pos(i,n+1),pos(i,n)) = -1;
%         Q(pos(i,n+1)) = 0;
    end
end


spy(F);
F = sparse(F);
T = F\Q;
T = reshape(T,n+2,n+2);

function ij = pos(i,j)
    N = size(v,1)+2;
    i=i+1; j=j+1; %omdat er een extra boord toegevoegd is
    ij = i+N.*(j-1);
end
end