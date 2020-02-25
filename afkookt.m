function cost_p = afkookt(v)

[T, K] = grid_discretisatie_khoek(v);

n = round(sqrt(length(v))); v = reshape(v,n,n); 
n = size(v,1) + 1;
p = 1;
k = v^p*65 + (1 - v^p)*0.2;
T = T(:);

lambda = K'\cost_T(T)';
T = reshape(T,n,n);

DK = zeros(n^2,(n-1)^2);

for i = 2:n-2
    for j = 2:n-2
        north = dmean(k(i,j),k(i-1,j));
        south = dmean(k(i,j),k(i+1,j));
        east = dmean(k(i,j),k(i,j+1));
        west = dmean(k(i,j),k(i,j-1));
        
        dK=[east+north, -east, 0, -north;
            -west, west+north, -north, 0;
            0, -south, south+west, -west;
            -south, 0, -east, south+east;
           ];
        
        Tij = [T(i+1,j) T(i+1,j+1) T(i,j+1) T(i,j)]';
        col = dK*Tij;
        
        DK(pos(i+1,j),posk(i,j)) = col(1);
        DK(pos(i+1,j+1),posk(i,j)) = col(2);
        DK(pos(i,j+1),posk(i,j)) = col(3);
        DK(pos(i,j),posk(i,j)) = col(4);
    end
end

cost_p = -lambda'*DK;
cost_p = cost_p(:);
% cost_p = reshape(cost_p,n-1,n-1);


    function c = cost_T(T)
        c = ((T.^15).*(sum(T.^16)).^(-15/16))';
    end
     
    function dk = dmean(k1, k2)
       dk = 2*k2^2/(k1+k2)^2; 
    end

    function ij = posk(i,j)
        N = size(v,1);
        ij = i+N.*(j-1);
    end

    function ij = pos(i,j)
        N = size(v,1) + 1;
        ij = i+N.*(j-1);
    end
end