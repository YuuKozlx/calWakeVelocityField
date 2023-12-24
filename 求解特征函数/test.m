% 打靶法求解

clc;
clear;
load '..\sigmoid函数产生NV频率\hp.mat'
load '..\sigmoid函数产生NV频率\Np.mat'
load ..\求解色散关系\'色散关系 ω_k 结果'\cp.mat

i = 1;
for guess = 0.001
    y_end(i, :) = shootingMethod(20, guess);
    i = i + 1;
end

% [x,y] = shootingMethod(0,0.001);

function y_end = shootingMethod(in_k, in_yp_guess)
% function [x,y] = shootingMethod(in_k,in_yp_guess)
k = in_k;
yp_guess = in_yp_guess;
cp = 0.0202303900200414;
% 定义变系数函数
N = @(x) Nv(x);
k3 = @(x) N(x)^2 / cp^2 - k^2;


% 定义微分方程
f = @(x, y, yp) yp; % y' = yp
g = @(x, y, yp) -k3(x) * y; % y'' = -b(x)y' - c(x)y

% 定义边值条件
y0 = 0; % 初始值
yN = 0; % 终止值

% 设置求解区域
a_val = 0;
b_val = 0.5;

% 设置步长和初始猜测
h = 0.05;
% yp_guess = 0.34;

% 打靶法求解
[x, y] = ode45(@(x, y) [f(x, y(1), y(2)); g(x, y(1), y(2))], [a_val, b_val], [y0; yp_guess]);

y_end = y(end, 1);

% 绘制结果
figure(1);
plot(y(:, 1),x, '-o');
xlabel('x');
ylabel('y(x)');
title('Shooting Method for BVP of Second-Order ODE with Variable Coefficients');
set(gca,'YDir','reverse');
end


function shootingMethodVariableCoeff()

% 定义变系数函数
a = @(x) x;
b = @(x) 1;
c = @(x) -x;

% 定义微分方程
f = @(x, y, yp) yp; % y' = yp
g = @(x, y, yp) -b(x) * yp - c(x) * y; % y'' = -b(x)y' - c(x)y

% 定义边值条件
y0 = 0; % 初始值
yN = 0; % 终止值

% 设置求解区域
a_val = 0;
b_val = 5;

% 设置步长和初始猜测
h = 0.1;
yp_guess = 1;

% 打靶法求解
[x, y] = ode45(@(x, y) [f(x, y(1), y(2)); g(x, y(1), y(2))], [a_val, b_val], [y0; yp_guess]);

% 绘制结果
figure(2);
plot(x, y(:, 1), '-o');
xlabel('x');
ylabel('y(x)');
title('Shooting Method for BVP of Second-Order ODE with Variable Coefficients');

end

function N = Nv(z)
N = zeros(size(z));

condition1 = (z >= 0 & z < 0.138);
N(condition1) = 0.010;
condition2 = (z >= 0.138 & z < 0.176);
N(condition2) = 0.400;
condition3 = (z >= 0.176 & z < 0.190);
N(condition3) = 0.600;
condition4 = (z >= 0.190 & z < 0.202);
N(condition4) = 0.800;
condition5 = (z >= 0.202 & z < 0.212);
N(condition5) = 1;
condition6 = (z >= 0.212 & z < 0.222);
N(condition6) = 1.200;
condition7 = (z >= 0.222 & z < 0.232);
N(condition7) = 1.400;
condition8 = (z >= 0.232 & z < 0.268);
N(condition8) = 1.600;
condition9 = (z >= 0.268 & z < 0.279);
N(condition9) = 1.400;
condition10 = (z >= 0.279 & z < 0.289);
N(condition10) = 1.200;
condition11 = (z >= 0.289 & z < 0.299);
N(condition11) = 1;
condition12 = (z >= 0.299 & z < 0.310);
N(condition12) = 0.800;
condition13 = (z >= 0.310 & z < 0.325);
N(condition13) = 0.600;
condition14 = (z >= 0.325 & z < 0.362);
N(condition14) = 0.400;
condition15 = (z >= 0.362 & z < 0.800);
N(condition15) = 0.010;

end

% 模态1 k = 0 打靶法 y'(0) = 0.0104