clc;
clear;
load '..\sigmoid函数产生NV频率\hp.mat'
load '..\sigmoid函数产生NV频率\Np.mat'
load ..\求解色散关系\'色散关系 ω_k 结果'\cp.mat
Np = Np';
hp = hp';

k_index = 61;
cp0 = cp(k_index, 2);
k = (k_index-1)*0.5;

% 定义边界条件
% syms alpha;
% syms beta;
alpha = 0;
beta = 0;
% 定义求解区域
a_val = 0;
b_val = 0.8;

% 设置网格参数
Nz = 100; % 网格点数
dz = (b_val - a_val) / (Nz - 1); % 网格步长
z = linspace(a_val, b_val, Nz);

% 初始化数值解数组
phi(1) = alpha;
phi(Nz) = beta;

% 通过二阶的有限差分法将微分方程化为关于网格点的线性方程组进行求解
% 通过递推求解系数矩阵
zz = z(2:end - 1);
q = dz ^ 2 * (Nv(zz) .^ 2 / cp0 ^ 2 - k ^ 2) - 2;

C(1, :) = [q(1), 1, zeros(1, Nz - 4)];

for i = 2:Nz - 2

    if i == Nz - 2
        C(i, Nz - 3:Nz - 2) = [1, q(i)];
    else
        C(i, i - 1:i + 1) = [1, q(i), 1];
    end

end

b = [alpha;zeros(Nz - 4, 1);beta];
det(C);

y = linsolve(C, b);
figure(1)
plot(y,zz);
set(gca, 'YDir', 'reverse');

% Nv频率
function N = Nv(z)
    N = zeros(size(z));

    condition1 = (z >= 0 & z < 0.138);
    N(condition1) = 0.0;
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
    N(condition15) = 0.0;

end
