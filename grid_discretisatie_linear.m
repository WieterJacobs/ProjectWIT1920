function [T, K] = grid_discretisatie_linear(v, p)
% K*T = Q
% v is verhouding metaal/vak (moet als dimensie een kwadraat hebben)
% p is exponent in SIMP methode

n = round(sqrt(length(v)));
v = reshape(v,n,n); 
n = size(v,1) + 1;

sizex = 0.01/2;
sizez = 0.001;
deltax = sizex/n;
k = v.^p*65 + (1 - v.^p)*0.2;
K = zeros((n)^2);
Q = zeros(n); Q = Q(:);

for i = 2:n-1
    for j = 2:n-1
        %Q(pos(i,j)) = 2*deltax^2/(2*sizex)^2/sizez;
        w = zeros(4,1);
        w(1) = mean(k(i-1,j-1),k(i-1,j));
        w(2) = mean(k(i,j-1),k(i,j));
        w(3) = mean(k(i-1,j-1),k(i,j-1));
        w(4) = mean(k(i-1,j),k(i,j));

        K(pos(i,j),pos(i-1,j)) = -w(1);
        K(pos(i,j),pos(i+1,j)) = -w(2);
        K(pos(i,j),pos(i,j-1)) = -w(3);
        K(pos(i,j),pos(i,j+1)) = -w(4);
        K(pos(i,j),pos(i,j)) = sum(w);
    end
end

%Dir = 1+round(.3*n):n-round(.3*n);
Dir = 1+round(0*n):n;

for j = 1:n
    K(pos(1,j),pos(1,j)) = 1;
    K(pos(1,j),pos(2,j)) = -1;
    K(pos(n,j),pos(n,j)) = 1;
    K(pos(n,j),pos(n-1,j)) = -1;
end
for i = 2:n-1
%     if ismember(i, Dir)
        K(pos(i,1),pos(i,1)) = 1;
        Q(pos(i,1)) = 0;
        K(pos(i,n),pos(i,n)) = 1;
        Q(pos(i,n)) = 10;
%     else
%         K(pos(i,1),pos(i,1)) = 1;
%         K(pos(i,1),pos(i,2)) = -1;
%end
    %K(pos(i,n),pos(i,n)) = 1;
    %K(pos(i,n),pos(i,n-1)) = -1;
end



% spy(F);
Ksparse = sparse(K);
T = Ksparse\Q;
% T = reshape(T,n,n);
% T = T(3:n-2,3:n-2);

    function k = mean(k1, k2)
        k = 2*k1*k2/(k1 + k2);
%         k = (k1 + k2)/2;
    end

    function ij = pos(i,j)
        N = size(v,1) + 1;
        ij = i+N.*(j-1);
    end
end