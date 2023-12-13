% 由于maltab无法直接实现多维矩阵按页乘法，因此使用mtimesx工具箱进行矩阵乘法
% mtimesx工具箱下载地址：https://ww2.mathworks.cn/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support
% mtimesx工具箱安装方法：将下载的文件夹放在matlab安装目录下的toolbox文件夹下，然后在matlab命令行中输入addpath(genpath('../mtimesx'))，即可完成安装
% mtimesx工具箱使用方法：在matlab命令行中输入help mtimesx，即可查看使用方法
% github地址：https://github.com/cybertk/mtimesx
addpath(genpath('../mtimesx'))

% 以下程序进行了矩阵化处理,对比soluteDispersionRelation.m程序，可以发现，后者程序针对每一个kp和omega都要进行单独进行判断
% 且针对 F = a0*a1*a2...*an 这种形式的矩阵乘法,计算ai的4个元素时首先进行判断,在进行串行计算.计算完成后,还要重组为矩阵
% 再进行矩阵乘法计算,耗时比较长
% 本程序则是将则是提前进行判断,将判断结果存储在矩阵中,然后进行多维矩阵乘法按页计算,耗时约为原来的1/2
% 事实上,也可以自行按照矩阵乘法规则编写矩阵乘法程序,将判断结果存储在矩阵中,然后进行矩阵乘法计算,这样耗时约为原来的1/10
% 但是,这种办法计算结果似乎和内置矩阵乘法和mtimesx工具箱计算结果有一定的差异,原因不明.若能完成优化,则可以进一步提高计算速度.
clc;
clear;
format long;
temp_load = load('../函数产生NV频率/hp.mat'); % 量化后的NV频率数据对应的深度，对应于论文中的hp
hp = temp_load.hp;
hp = hp';
temp_load = load('../函数产生NV频率/Np.mat'); % 量化后的NV频率数据对应的Np，对应于论文中的Np
Np = temp_load.Np;
NNp = Np(1:end-1);
NNp = NNp';
w = 0.00001:0.0002:1.60001;
w_repmat = repmat(w, length(NNp), 1);

mpIsRe = (w_repmat > NNp);
mpIsIm = 1 - mpIsRe;

f12 = zeros(60, length(w));


tStart = tic;

kpindex = 1;
max_kp = 30;

for kp = 0.5:0.5:max_kp

    omega = w;

    k3p = sqrt(kp .^ 2 .* (NNp .^ 2 ./ omega .^ 2 - 1));
    mp = k3p * 1j;

    F_temp(:, :, 1, 1) = mpIsRe .* cosh(mp .* hp) + mpIsIm .* cos(k3p .* hp);

    F_temp(:, :, 1, 2) = mpIsRe .* sinh(mp .* hp) ./ mp + mpIsIm .* sin(k3p .* hp) ./ k3p;

    F_temp(:, :, 2, 1) = mpIsRe .* mp .* sinh(mp .* hp) - mpIsIm .* sin(k3p .* hp) .* k3p;

    F_temp(:, :, 2, 2) = mpIsRe .* cosh(mp .* hp) + mpIsIm .* cos(k3p .* hp);

    F_temp_permuted = permute(F_temp, [3, 4, 1, 2]);
    F = F_temp_permuted(:, :, 1, :);

    for j = 2:length(NNp)
        F = mtimesx(F_temp_permuted(:, :, j, :), F);
    end

    f12(kpindex, :) = (squeeze(F(1, 2, :)))';

    if mod(kpindex, 1000) == 0 && kp ~= max_kp
        F12 = f12(kpindex - 999:kpindex, :);
        save(['f12_', num2str(kpindex - 999), '_', num2str(kpindex), '.mat'], 'F12');
    elseif kp == max_kp
        F12 = f12(kpindex - mod(kpindex, 1000) + 1:kpindex, :);
        save(['f12_', num2str(kpindex - mod(kpindex, 1000) + 1), '_', num2str(kpindex), '.mat'], 'F12');
    end

    kpindex = kpindex + 1;
end

clear a11_mpIsRe a11_mpIsIm a12_mpIsRe a12_mpIsIm a21_mpIsRe a21_mpIsIm a22_mpIsRe a22_mpIsIm apJudge
toc(tStart)

