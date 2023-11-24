clc;
clear;
format long;

load('data_1_SaltWaterDensity.mat');

[density, depth] = saltWaterDensityInterpolation(data_1_SaltWaterDensity);
figure(1);
plot(data_1_SaltWaterDensity(:, 1), data_1_SaltWaterDensity(:, 2), 'o', depth, density, '-');
title('盐水密度数据插值');
xlabel('盐水密度');
ylabel('深度');
set(gca, 'YDir', 'reverse');






% 三次样条给盐水密度数据插值，得到盐水分层数据，包括盐水密度和深度数据
function [density, depth] = saltWaterDensityInterpolation(originalData)
    x = originalData(:, 2); % x 坐标是 密度数据
    y = originalData(:, 1); % y 坐标是 深度数据
    pp = spline(x, y);
    xx = min(x):0.001:max(x);
    yy = ppval(pp, xx);
    depth = yy;
    density = xx;
end

% 计算盐水密度分布的NV频率
function [nv] = nvFrequency(density, depth)
    g = 9.8;
    rho0 = 1000;
    rho = density;
    z = depth;
    nv = sqrt(g * (rho0 - rho) ./ rho .* z);
end
