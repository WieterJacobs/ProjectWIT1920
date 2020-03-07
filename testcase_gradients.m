v = zeros(20); v = v(:);
p = 1;
cost1 = afkookt(v,p);
cost1 = reshape(cost1,20,20);
figure; heatmap(cost1);
cost2 = fin_dif(v(:),p);
cost2 = reshape(cost2,20,20);
figure; heatmap(cost2);
figure; heatmap((cost2-cost1)./(cost2+1));
