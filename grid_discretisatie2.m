function T = grid_discretisatie2(v)
% F*T = Q
n = size(v,1);
vnonce = v;
v = zeros(n+4);
v(3:n+2,3:n+2) = vnonce;
n = size(v,1);

deltax = 0.01;
deltan = deltax/n;
k = v*65 + (1 - v)*0.2;
F = zeros((n)^2);
Q = zeros(n); Q = Q(:);

for i = 2:n-1
    for j = 2:n-1
        Q(pos(i,j)) = 2*deltan^2*10^8; % !!!
        w = zeros(4,1);
        w(1) = (k(i-1,j)+k(i,j))/2;
        w(2) = (k(i+1,j)+k(i,j))/2;
        w(3) = (k(i,j-1)+k(i,j))/2;
        w(4) = (k(i,j+1)+k(i,j))/2;

        F(pos(i,j),pos(i-1,j)) = -w(1);
        F(pos(i,j),pos(i+1,j)) = -w(2);
        F(pos(i,j),pos(i,j-1)) = -w(3);
        F(pos(i,j),pos(i,j+1)) = -w(4);
        F(pos(i,j),pos(i,j)) = sum(w);
    end
end

Dir = 1+round(.3*n):n+2-round(.3*n);

for j = 1:n
    F(pos(1,j),pos(1,j)) = 1;
    F(pos(1,j),pos(2,j)) = -1;
    F(pos(n,j),pos(n,j)) = 1;
    F(pos(n,j),pos(n-1,j)) = -1;
end
for i = 2:n-1
    if ismember(i, Dir)
        F(pos(i,1),pos(i,1)) = 1;
        Q(pos(i,1)) = 293;
        F(pos(i,n),pos(i,n)) = 1;
        Q(pos(i,n)) = 293;
    else
        F(pos(i,1),pos(i,1)) = 1;
        F(pos(i,1),pos(i,2)) = -1;
        F(pos(i,n),pos(i,n)) = 1;
        F(pos(i,n),pos(i,n-1)) = -1;
    end
end


spy(F);
F = sparse(F);
T = F\Q;
T = reshape(T,n,n);
% T = T(3:n-2,3:n-2);

function ij = pos(i,j)
    N = size(v,1);
    ij = i+N.*(j-1);
end
end