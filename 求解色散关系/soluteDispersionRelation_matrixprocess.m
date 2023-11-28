% 进行矩阵化处理 速度快了两倍

addpath(genpath('../mtimesx'))
clc;
clear;
format long;
temp_load = load('../sigmoid函数产生NV频率/hp.mat'); % 量化后的NV频率数据对应的深度，对应于论文中的hp
hp = temp_load.hp;
hp = hp';
temp_load = load('../sigmoid函数产生NV频率/Np.mat'); % 量化后的NV频率数据对应的Np，对应于论文中的Np
Np = temp_load.Np;
NNp = Np(1:15);
NNp = NNp';
w = 0.00001:0.00001:1.6;
w_repmat = repmat(w, length(NNp), 1);

mpIsRe = (w_repmat > NNp);
mpIsIm = 1 - mpIsRe;

f12 = zeros(20, length(w));

tStart = tic;

kpindex = 1;
max_kp = 0.51;
for kp = 0.005:0.005:max_kp

    omega = 0.000005:0.00001:1.6;

    k3p = sqrt(kp.^2.*(NNp.^2 ./ omega.^2 - 1));
    mp = k3p * 1j;

    F_temp(:, :, 1, 1) = mpIsRe .* cosh(mp.*hp) + mpIsIm .* cos(k3p.*hp);

    F_temp(:, :, 1, 2) = mpIsRe .* sinh(mp.*hp) ./ mp + mpIsIm .* sin(k3p.*hp) ./ k3p;

    F_temp(:, :, 2, 1) = mpIsRe .* mp .* sinh(mp.*hp) - mpIsIm .* sin(k3p.*hp) .* k3p;

    F_temp(:, :, 2, 2) = mpIsRe .* cosh(mp.*hp) + mpIsIm .* cos(k3p.*hp);

    F_temp_permuted = permute(F_temp, [3, 4, 1, 2]);
    F = F_temp_permuted(:, :, 1, :);

    for j = 2:15
        F = mtimesx(F_temp_permuted(:, :, j, :), F);
    end

    f12(kpindex, :) = (squeeze(F(1, 2, :)))';

    if mod(kpindex, 1000) == 0 && kp ~= max_kp
        F12 = f12(kpindex-999:kpindex, :);
        save(['f12_', num2str(kpindex-999), '_', num2str(kpindex), '.mat'], 'F12');
    elseif kp == max_kp
        F12 = f12(kpindex-mod(kpindex, 1000)+1:kpindex, :);
        save(['f12_', num2str(kpindex-mod(kpindex, 1000)+1), '_', num2str(kpindex), '.mat'], 'F12');
    end

    kpindex = kpindex + 1;
end

clear a11_mpIsRe a11_mpIsIm a12_mpIsRe a12_mpIsIm a21_mpIsRe a21_mpIsIm a22_mpIsRe a22_mpIsIm apJudge
toc(tStart)


