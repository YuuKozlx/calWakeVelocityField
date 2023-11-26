
w_index = find_sign_change(f12);
w = 0.00001*w_index;
k = 5:5:100;
w = fliplr(w);
figure(1)
p = plot(k,w);
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
xlabel('k','FontSize',14,'FontName','Times New Roman');
xlim([0, 110])
xticks(0:20:120) % 设置 X 轴的大刻度分度值，范围是 0 到 110，步长为 10
ylim([0, 1.8])
yticks(0:0.2:1.8) % 设置 Y 轴的大刻度分度值，范围是 0 到 1.8，步长为 0.2
ylabel('\omega', 'FontSize', 14, 'FontName', 'Times New Roman');
title('The relationship between \omega and k','FontSize',14,'FontName','Times New Roman');
legend('model1','model2','model3','Location','southeast');
grid on;
set(gca,'FontSize',14,'XMinorTick','on','YMinorTick','on');
set(gcf,'Position',[100 100 500 400]);

function result_indices = find_sign_change(A)

    % 初始化异号点的位置矩阵
    result_indices = [];

    % 找出前后两点异号的位置或等于0的位置
    for i = 1:size(A, 1)
        row = A(i, :);
        sign_changes = find((row(1:end - 1) .* row(2:end) < 0)); % 修改这里
        zero_indices = find(abs(row) <= 1e-15);

        result = [sign_changes, zero_indices];
        result = sort(result);

        if length(result) >= 4 % 只取前4个
            result = result(end-3:end);
        end

        result_indices = [result_indices; result];
    end

end
