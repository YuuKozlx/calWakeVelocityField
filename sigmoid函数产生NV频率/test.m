% 定义连续信号
t = linspace(0, 1, 1000);
x = 3*sin(2 * pi * 5 * t);

% 设置多级阈值
thresholds = [0.2, 0.5, 1.0];

% 设置初始步长
initial_step_size = 1;

% 进行多级自适应量化
y_quantized = adaptiveQuantizationMultiLevel(x, thresholds, initial_step_size);

% 绘制原始信号和多级自适应量化后的信号
figure;
subplot(2, 1, 1);
plot(t, x, 'b', 'LineWidth', 2);
title('原始信号');
xlabel('时间');
ylabel('幅度');

subplot(2, 1, 2);
stairs(t, y_quantized, 'r', 'LineWidth', 2);
title('多级自适应量化后的信号');
xlabel('时间');
ylabel('量化值');



function y_quantized = adaptiveQuantizationMultiLevel(x, thresholds, initial_step_size)
    % x: 输入信号
    % thresholds: 幅值变化的多级阈值
    % initial_step_size: 初始步长

    % 初始化
    y_quantized = zeros(size(x));
    step_size = initial_step_size;  % 初始步长

    % 进行多级自适应量化
    for i = 2:length(x)
        % 计算当前点与前一个点的差值
        delta = x(i) - x(i-1);

        % 根据差值选择相应的阈值和步长
        for j = 1:length(thresholds)
            if abs(delta) > thresholds(j)
                step_size = initial_step_size / (2^j);  % 调整步长
            end
        end

        % 进行量化
        y_quantized(i) = round(x(i) / step_size) * step_size;
    end
end

