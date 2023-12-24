clc;
clear;
format long;

% 深度序列
depth = 0:0.001:0.8;
% 采用sigmoid函数产生密度序列
density = densityDistribution(depth);
figure(1);
plot(density, depth, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.5);
title('盐水密度分布');
xlabel('密度');
ylabel('深度');
set(gca, 'YDir', 'reverse');

figure(2);
nvFreq = nvFrequency(depth, density);
plot(nvFreq, depth, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.5);
title('盐水密度分布的NV频率');
xlabel('NV频率');
ylabel('深度');
set(gca, 'YDir', 'reverse');

figure(3)
% 绘制量化后的频率序列
nvQuant = nvQuantization(nvFreq);
plot(nvQuant, depth, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.5);
title('盐水密度分布的NV频率');
xlabel('NV频率');
ylabel('深度');
set(gca, 'YDir', 'reverse', 'YLim', [0, 0.8], 'XLim', [0, 1.8]);
save('nvQuant.mat', 'nvQuant');

figure(4)
% 绘制量化后频率发生变化的深度和对应的NV频率
depthQuantIndex = nvDepth(nvQuant);
depthMark = depth(depthQuantIndex);
hp = diff(depthMark);
Np = nvQuant(depthQuantIndex);
save('hp.mat', 'hp');
% writematrix(hp, 'hp.csv');
save('Np.mat', 'Np');
stem(depth(depthQuantIndex), nvQuant(depthQuantIndex), 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.5);
xlabel('深度');
ylabel('NV频率');

% 产生盐水密度分布
function [density] = densityDistribution(depth)
% 使用sigmoid函数产生密度分布
% sigmoid(x) = 1 / (1 + exp(-x))
k = 40;
z0 = 0.3;
density = 22.5 ./ (1 + exp(-k*(depth - z0))) + 1000;
end

% 计算盐水密度分布的NV频率和深度的关系
% 浮力频率计算公式：nv = sqrt(g * rho0'/ rho0);
function [nv] = nvFrequency(depth, density)
% 使用sigmoid函数计算NV频率
% sigmoid函数的导数为
% sigmoid'(x) = sigmoid(x) * (1 - sigmoid(x))
z = depth;
% 密度分布的导数
% densityDerivative = (1125 * exp(-10 * (z - 0.25))) ./ (exp(-10 * (z - 0.25)) + 1) .^ 2;
k = 40;
z0 = 0.3;
densityDerivative = 22.5 * k .* exp(-k.*(z - z0)) ./ (1 + exp(-k.*(z - z0))).^2;

nv = sqrt(9.8*densityDerivative./density);
end

% NV频率的自适应量化
function [nvQuant] = nvQuantization(nvFreq)
% 根据多级阈值进行自适应量化
% 深度步长
initial_step_size = 0.1;
% 量化的阈值
threshold = [0.005, 0.01, 0.02, 0.03, 0.04]*(max(nvFreq)-min(nvFreq));
nvQuant = zeros(size(nvFreq));
% 进行多级自适应量化
for i = 2:length(nvFreq)
    % 计算当前深度的量化值
    delta = abs(nvFreq(i)-nvFreq(i - 1));

    % 根据插值选择量化阈值和步长
    % 判断当前深度的NV频率是否超过阈值
    for j = 1:length(threshold)

        if delta < threshold(j)
            % 选择当前阈值和步长
            step_size = initial_step_size / (2^(j - 1));
            break;
        end

    end

    % 量化当前深度的NV频率
    nvQuant(i) = round(nvFreq(i)/step_size) * step_size;
end
    nvQuant(nvQuant==0) = 0.01;

end

% 在NV频率量化后，记录NV频率恰发生变化的深度(包括起始深度和结束深度)，并返回索引。
function [depthQuantIndex] = nvDepth(nvQuant)
% 记录NV频率发生变化的深度
depthQuantIndex = 1;

for i = 2:length(nvQuant)

    if nvQuant(i) ~= nvQuant(i - 1)
        depthQuantIndex = [depthQuantIndex, i];
    end

end

depthQuantIndex = [depthQuantIndex, length(nvQuant)];

end
