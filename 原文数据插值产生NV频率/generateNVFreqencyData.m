clc;
clear;
format long;

load('data_1_SaltWaterDensity.mat');

[density, depth] = saltWaterDensityInterpolation(data_1_SaltWaterDensity);

ans = diff(density);

figure(1);
plot(data_1_SaltWaterDensity(:, 1), data_1_SaltWaterDensity(:, 2), 'o', density, depth, '-');
title('盐水密度数据插值');
xlabel('盐水密度');
ylabel('深度');
set(gca, 'YDir', 'reverse');

nv = nvFrequency(density, depth);
figure(2);
plot(nv, depth, '-');
title('盐水密度分布的NV频率');
xlabel('NV频率');
ylabel('深度');
set(gca, 'YDir', 'reverse');

% 三次样条给盐水密度数据插值，得到盐水分层数据，包括盐水密度和深度数据
function [density, depth] = saltWaterDensityInterpolation(originalData)
    x = originalData(:, 2); % x 坐标是 深度数据
    y = originalData(:, 1); % y 坐标是 密度数据
    pp = spline(x, y);
    xx = 0:0.001:0.8;
    yy = ppval(pp, xx);
    density = yy;
    depth = xx;
    % 保留两位小数
    density = round(density, 4);
    depth = round(depth, 4);
end

% 计算盐水密度分布的NV频率
% 浮力频率计算公式：nv = sqrt(g * rho0'/ rho0);
function [nv] = nvFrequency(density, depth)
    g = 9.8;
    rho = density;
    z = depth;
    diffRho = diff(rho);
    diffZ = diff(z);
    rho_Prime = diffRho ./ diffZ;
    rho_Prime = [0, rho_Prime];
    nv = sqrt(g * rho_Prime ./ rho);
end
