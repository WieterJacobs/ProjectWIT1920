function T = grid_discretisatie(v)
% F*T = Q
n = size(v,1);
vnonce = v;
v = zeros(n+4);
v(3:n+2,3:n+2) = vnonce;
n = size(v,1);

deltax = 0.01;
deltan = deltax/n;
k = v*65 + (1 - v)*0.2;
F = zeros(n^2);
Q = zeros(n); Q = Q(:);

for i = 3:n-2
    for j = 3:n-2
        Q(pos(i,j)) = 2*deltan^2 *10^7; %!!!
        
        F(pos(i,j),pos(i,j)) = 16*k(i,j) + k(i-1,j) + k(i+1,j) +  k(i,j-1) + k(i,j+1);
        F(pos(i,j),pos(i-1,j)) = -4*k(i,j) -4*k(i-1,j);
        F(pos(i,j),pos(i+1,j)) = -4*k(i,j) -4*k(i+1,j);
        F(pos(i,j),pos(i,j-1)) = -4*k(i,j) -4*k(i,j-1);
        F(pos(i,j),pos(i,j+1)) = -4*k(i,j) -4*k(i,j+1);
        
        F(pos(i,j),pos(i-1,j-1)) = k(i-1,j) + k(i,j-1);
        F(pos(i,j),pos(i-1,j+1)) = k(i-1,j) + k(i,j+1);       
        F(pos(i,j),pos(i+1,j-1)) = k(i+1,j) + k(i,j-1);
        F(pos(i,j),pos(i+1,j+1)) = k(i+1,j) + k(i,j+1);
        
        F(pos(i,j),pos(i-2,j)) = k(i-1,j);
        F(pos(i,j),pos(i+2,j)) = k(i+1,j);
        F(pos(i,j),pos(i,j-2)) = k(i,j-1);
        F(pos(i,j),pos(i,j+2)) = k(i,j+1);
    end
end

Dir = 1+round(.3*n):n+2-round(.3*n);
for j = 1:n
    for i = [1 2]
        F(pos(i,j),pos(i,j)) = 1;
        F(pos(i,j),pos(i+1,j)) = -1;
        F(pos(n+1-i,j),pos(n+1-i,j)) = 1;
        F(pos(n+1-i,j),pos(n-i,j)) = -1;
    end
end
for i = 2:n-1
    for j = [1 2]
        if ismember(i, Dir)
            F(pos(i,j),pos(i,j)) = 1;
            Q(pos(i,j)) = 293;
            F(pos(i,n+1-j),pos(i,n+1-j)) = 1;
            Q(pos(i,n+1-j)) = 293;
        else
            F(pos(i,j),pos(i,j)) = 1;
            F(pos(i,j),pos(i,j+1)) = -1;
            F(pos(i,n+1-j),pos(i,n+1-j)) = 1;
            F(pos(i,n+1-j),pos(i,n-j)) = -1;
        end
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