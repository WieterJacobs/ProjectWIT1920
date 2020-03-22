function cost_p = afkookt(v,p)

v = v(:);
kder = p*v.^(p-1)*65 + (1 - p*v.^(p-1))*0.2;
[T, K] = grid_discretisatie_khoek(v, p);

n = round(sqrt(length(v))); v = reshape(v,n,n); 
n = size(v,1) + 1;
k = v^p*65 + (1 - v^p)*0.2;
T = T(:);

lambda = K'\cost_T(T)';
T = reshape(T,n,n);

DK = zeros(n^2,(n-1)^2);

for i = 1:n-1
    for j = 1:n-1
        try
            north = dmean(k(i,j),k(i-1,j));
        catch
            north = 0;
        end
        try
            south = dmean(k(i,j),k(i+1,j));
        catch
            south = 0;
        end
        try
            east = dmean(k(i,j),k(i,j+1));
        catch
            east = 0;
        end
        try
            west = dmean(k(i,j),k(i,j-1));
        catch
            west = 0;
        end
            
      
%         dK=[east+north, -east, 0, -north;
%             -west, west+north, -north, 0;
%             0, -south, south+west, -west;
%             -south, 0, -east, south+east;
%            ];
        dK=[south+west, -south, 0, -west;
            -south, east+south, -east, 0;
            0, -east, east+north, -north;
            -west, 0, -north, west+north;
           ];
       
        Tij = [T(i+1,j) T(i+1,j+1) T(i,j+1) T(i,j)]';
        col = dK*Tij;
        
        Dir = 1+round(.6*n):n+2;
        
        if ~(j == 1 && ismember(i+2,Dir))
            DK(pos(i+1,j),posk(i,j)) = col(1);
        end
        DK(pos(i+1,j+1),posk(i,j)) = col(2);
        DK(pos(i,j+1),posk(i,j)) = col(3);
        if ~(j == 1 && ismember(i+1,Dir))
            DK(pos(i,j),posk(i,j)) = col(4);
        end
    end
end

cost_p = -lambda'*DK;
cost_p = cost_p(:).*kder;
% cost_p = reshape(cost_p,n-1,n-1);


    function c = cost_T(T)
        c = ((T.^15).*(sum(T.^16)).^(-15/16))';
    end
     
    function dk = dmean(k1, k2)
       dk = 2*k2^2/(k1+k2)^2; 
%        dk = .5; 
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