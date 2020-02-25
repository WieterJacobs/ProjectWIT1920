v = zeros(20);
cost1 = afkookt(v(:));
figure; heatmap(cost1);
cost2 = fin_dif(v);
figure; heatmap(cost2);