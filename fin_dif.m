function cost_p = fin_dif(v,p)
v = reshape(v,sqrt(length(v)),sqrt(length(v)));
cost_p = zeros(size(v,1)); %
for i = 1:length(v(:))
    h = 0.0001;
    vv = v;
    
    vv(i) = v(i) + h;
    T1 = grid_discretisatie_khoek(vv(:),p);
    vv(i) = v(i) - h;
    T2 = grid_discretisatie_khoek(vv(:),p);
    cost_p(i) = (cost_T(T1) - cost_T(T2))/2/h;
end
cost_p = cost_p(:);
end

    function c = cost_T(T)
        c = sum(T(:).^16)^(1/16);
    end
