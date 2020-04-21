function [error,rico] = discretisation_simulation(m,p,method) %m is number 
                                                        %of grids to be tested, 
                                                        %method is to test
                                                        %either sin, sin2
                                                        %or sincos
    gridsize = linspace(21,20*(m+1)+1,m);
    result = zeros(m,1);
    ricos = zeros(m-1,1);
    index = 1;
    testx = rand;
    testy = rand;
    if method == "sin"
        for n = gridsize
           v = ones(n);
           T = reshape(grid_discretisatie_sin(v(:),p),n+1,n+1);
           for i = 1:n+1
               for j = 1:n+1
                   T(i,j) = abs(T(i,j)-(2+sin(2*pi*(j-0.5)/(n+1))))/(2+sin(2*pi*(j-0.5)/(n+1)));
               end
           end
           result(index) = T(1+floor(testx*n),1+floor(testy*n));
           index = index + 1;
        end
    elseif method == "sin2"
        for n = gridsize
           v = ones(n);
           T = reshape(grid_discretisatie_sin2(v(:),p),n+1,n+1);
           for i = 1:n+1
               for j = 1:n+1
                   T(i,j) = abs(T(i,j)-(2+sin(2*pi*(j-0.5)/(n+1))))/(2+sin(2*pi*(j-0.5)/(n+1)));
               end
           end
           result(index) = T(1+floor(testx*n),1+floor(testy*n));
           index = index + 1;
        end
    elseif method == "linear"
        for n = gridsize
           v = ones(n);
           T = reshape(grid_discretisatie_linear(v(:),p),n+1,n+1);
           for i = 1:n+1
               for j = 1:n+1
                   T(i,j) = abs(T(i,j)-10*(j-0.5)/(n+1))/10*(j-0.5)/(n+1);
               end
           end
           result(index) = T(1+floor(testx*n),1+floor(testy*n));
           index = index + 1;
        end
    end
    
    %loglog(gridsize,result);
    for i = 1:m-1
        tel = log(result(i+1))-log(result(i));
        noe = log(gridsize(i+1))-log(gridsize(i));
        ricos(i) = tel/noe;
    end
    error = result;
    rico = ricos;
end