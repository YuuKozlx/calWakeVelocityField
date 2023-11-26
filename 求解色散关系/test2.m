% 假设 x 和 y 是你的数据
w = (1:1:length(f12))*0.00001;
plot(w(10000:end), f12(10000:end), 'bo');
hold on; % 保持当前图像

% 绘制 y=0 的基准线
line([min(w) max(w)], [0 0], 'Color', 'red', 'LineStyle', '--');

hold off; % 关闭 hold on
