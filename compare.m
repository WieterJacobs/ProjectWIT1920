function [T_ok,K_ok,DK_ok,cost_p_ok] = compare(path_T,path_K,path_DK,path_cost_p)
    K_1=readmatrix(path_K);
    T_1=readmatrix(path_T);
    DK_1=readmatrix(path_DK);
    cost_p_1=path_cost_p;
    v=0.4*ones(4);
    [T_2, K_2] = grid_discretisatie_khoek(v(:));
    [cost_p_2,DK_2]=afkookt(v(:));
    
    T_ok=max(T_1-T_2)<10^(-3);
    K_ok=max(K_1-K_2)<10^(-3);
    DK_ok=max(DK_1-DK_2)<10^(-3);
    cost_p_ok=max(cost_p_1-cost_p_2)<10^(-3);
end
