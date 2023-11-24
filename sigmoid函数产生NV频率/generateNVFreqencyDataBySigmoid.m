clc;
clear;
format long;

% 深度序列
depth = 0:0.001:0.8;
% 采用sigmoid函数产生密度序列
density = densityDistribution(depth);
figure(1);
plot(density, depth, '-');
title('盐水密度分布');
xlabel('密度');
ylabel('深度');
set(gca, 'YDir', 'reverse');

figure(2);
nvFreq = nvFrequency(depth, density);
plot(nvFreq, depth, '-');
title('盐水密度分布的NV频率');
xlabel('NV频率');
ylabel('深度');
set(gca, 'YDir', 'reverse');

figure(3)
% 绘制量化后的频率序列
nvQuant = nvQuantization(depth);
plot(nvQuant, depth, '-');
title('盐水密度分布的NV频率');
xlabel('NV频率');
ylabel('深度');
set(gca, 'YDir', 'reverse');

% 产生盐水密度分布
function [density] = densityDistribution(depth)
    % 使用sigmoid函数产生密度分布
    % sigmoid(x) = 1 / (1 + exp(-x))
    density = 22.5 ./ (1 + exp(-50 * (depth - 0.25))) + 1000;
end

% 计算盐水密度分布的NV频率和深度的关系
% 浮力频率计算公式：nv = sqrt(g * rho0'/ rho0);
function [nv] = nvFrequency(depth, density)
    % 使用sigmoid函数计算NV频率
    % sigmoid函数的导数为
    % sigmoid'(x) = sigmoid(x) * (1 - sigmoid(x))
    z = depth;
    % 密度分布的导数
    densityDerivative = (1125 * exp(-50 * (z - 0.25))) ./ (exp(-50 * (z - 0.25)) + 1) .^ 2;

    nv = sqrt(9.8 * densityDerivative ./ density);
end

% 对计算得出的NV频率进行自适应量化，斜率越大，量化精度越高，斜率越小，量化精度越低
% function [nvQuant] = nvQuantization(depth)
%     z = depth;
%     nvFrequencyDerivative = (112500 * exp(25 - 100 * z)) ./ (exp(25/2 - 50 * z) + 1) .^ 3 - (56250 * exp(25/2 - 50 * z)) ./ (exp(25/2 - 50 * z) + 1) .^ 2;
%     nvFrequencyDerivative = abs(nvFrequencyDerivative);
%     nvFrequencyDerivative = nvFrequencyDerivative / max(nvFrequencyDerivative);
%     nvQuant = nvFrequencyDerivative;

% end

function [nvQuant] = nvQuantization(depth, nvFreq)
    % 根据多级阈值进行自适应量化
    % 量化步长
    step = 0.1;
    % 量化的阈值
    threshold = [0.1, 0.2, 0.3, 0.4, 0.5]*22.5;
    % 量化的级数
    level = length(threshold) + 1;
    
end
