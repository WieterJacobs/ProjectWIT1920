function [T_ok,K_ok,DK_ok,cost_p_ok] = compare(path_T,path_K,path_DK,path_cost_p)
    K_1=readmatrix(path_K);
    T_1=readmatrix(path_T);
    DK_1=readmatrix(path_DK);
    cost_p_1=readmatrix(path_cost_p);
    v=0.4*ones(4);
    [T_2, K_2] = grid_discretisatie_khoek(v(:));
    [cost_p_2,DK_2]=afkookt(v(:));
    %n = round(sqrt(length(K_1)));
    %cost_p_2 = reshape(cost_p_2,n-1,n-1);
    
    T_ok=max(T_1-T_2);
    K_ok=max(K_1-K_2);
    DK_ok=max(DK_1-DK_2);
    cost_p_ok=max(cost_p_1-cost_p_2);
end
