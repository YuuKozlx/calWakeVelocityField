w_index = find_sign_change(f12);
w = 0.00001 * w_index;
k = 5:5:100;
w = fliplr(w);
figure(1)
p = plot(k, w(:,1:3));
% p(4).LineWidth = 1.5;
% p(4).Marker = 'o';
% p(4).MarkerSize = 10;
% p(4).LineStyle = '-';
p(3).LineWidth = 1.5;
p(3).Marker = 'square';
p(3).MarkerSize = 10;
p(3).LineStyle = '-';
p(2).LineWidth = 1.5;
p(2).Marker = 'diamond';
p(2).MarkerSize = 10;
p(2).LineStyle = '-';
p(1).LineWidth = 1.5;
p(1).Marker = 'pentagram';
p(1).MarkerSize = 10;
p(1).LineStyle = '-';
xlabel('k', 'FontSize', 14, 'FontName', 'Times New Roman');
xlim([0, 110])
xticks(0:20:120) % 设置 X 轴的大刻度分度值，范围是 0 到 110，步长为 10
ylim([0, 1.8])
yticks(0:0.2:1.8) % 设置 Y 轴的大刻度分度值，范围是 0 到 1.8，步长为 0.2
ylabel('\omega', 'FontSize', 14, 'FontName', 'Times New Roman');
title('The relationship between \omega and k', 'FontSize', 14, 'FontName', 'Times New Roman');
legend('model1', 'model2', 'model3', 'Location', 'southeast');
grid on;
set(gca, 'FontSize', 14, 'XMinorTick', 'on', 'YMinorTick', 'on');
set(gcf, 'Position', [100 100 500 400]);





