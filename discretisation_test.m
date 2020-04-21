function [Rsin,Rsin2,Rlinear] = discretisation_test(m,p,n)
    Rsin = 0;
    for i = 1:n
        [T,R]=discretisation_simulation(m,p,"sin");
        Rsin = Rsin + mean(R);
    end
    Rsin = Rsin/n;
    Rsin2 = 0;
    for i = 1:n
        [T,R]=discretisation_simulation(m,p,"sin2");
        Rsin2 = Rsin2 + mean(R);
    end
    Rsin2 = Rsin2/n;
    Rlinear = 0;
    for i = 1:n
        [T,R]=discretisation_simulation(m,p,"linear");
        Rlinear = Rlinear + mean(R);
    end
    Rlinear = Rlinear/n;
end