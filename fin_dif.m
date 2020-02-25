function cost_p = fin_dif(v)
cost_p = zeros(size(v,1));
for i = 1:length(v(:))
    h = 0.0001;
    vv = v;
    
    vv(i) = v(i) + h;
    T1 = grid_discretisatie_khoek(vv);
    vv(i) = v(i) - h;
    T2 = grid_discretisatie_khoek(vv);
    cost_p(i) = cost_p(i) + (max(T1(:)) - max(T2(:)))/2/h;
end
end