function [T, K] = grid_discretisatie_sin2(v, p) %oplossing T=2+sin(2pi*x/deltax)
% K*T = Q
% v is verhouding metaal/vak (moet als dimensie een kwadraat hebben)
% p is exponent in SIMP methode

n = round(sqrt(length(v)));
v = reshape(v,n,n); 
k = zeros(n);
sizex = 0.01;
sizez = 0.001;
for i = 1:n
    for j = 1:n
        k(i,j) = ((i-0.5)*(j-0.5))*sizex^2/n^2;
    end
end
n = size(v,1) + 1;


deltax = sizex/n;
%k = v.^p*65 + (1 - v.^p)*0.2;
K = zeros((n)^2);
Q = zeros(n); Q = Q(:);

for i = 2:n-1
    for j = 2:n-1
        Q(pos(i,j)) = ((2*pi*(i-0.5)*deltax^2)/n)*(((2*pi*(j-0.5))/n)*sin((2*pi*(j-0.5))/n)-cos((2*pi*(j-0.5))/n));
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

Dir = 1:n;

for j = 1:n
    K(pos(1,j),pos(1,j)) = 1;
    K(pos(1,j),pos(2,j)) = -1;
    K(pos(n,j),pos(n,j)) = 1;
    K(pos(n,j),pos(n-1,j)) = -1;
end
for i = 2:n-1
%     if ismember(i, Dir)
        K(pos(i,1),pos(i,1)) = 1;
        Q(pos(i,1)) = 2;
        K(pos(i,n),pos(i,n)) = 1;
        Q(pos(i,n)) = 2;
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