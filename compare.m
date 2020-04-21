[T_diff,K_diff,cost_v_diff]=compare_matrices("v_matrix.txt","T_matrix.txt","K_matrix.txt","cost_v_matrix.txt",3);
T_diff
K_diff
heatmap(cost_v_diff);

function [T_ok,K_ok,cost_p_ok] = compare_matrices(path_v,path_T,path_K,path_cost_p,p)
    K_1=readmatrix(path_K);
    T_1=readmatrix(path_T);
    %DK_1=readmatrix(path_DK);
    cost_p_1=readmatrix(path_cost_p);
    v=readmatrix(path_v);
    [T_2, K_2] = grid_discretisatie_khoek(v(:),p);
    cost_p_2=afkookt(v(:),p);
    %n = round(sqrt(length(K_1)));
    cost_p_2=reshape(cost_p_2,20,20);
    T_ok=max(T_1-T_2);
    K_ok=max(max(K_1-K_2));
    %DK_ok=max(DK_1-DK_2);
    cost_p_ok=abs((cost_p_1-cost_p_2));
end
