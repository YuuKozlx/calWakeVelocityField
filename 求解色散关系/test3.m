clc; clear;
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
[row_w_eq_Np, col_w_eq_Np] = find(w_repmat == NNp);
apJudge = w_repmat > NNp;
clear w_repmat w
% apJudge(apJudge == 0) = 1;
%  w>Np, apJudge = 1 mp为实数
%  w<Np, apJudge = -1 mp为虚数
a11 = apJudge;
a12 = apJudge;
a21 = apJudge;
a22 = apJudge;
tic

for kp = 5:5:100
    % kp = 5;
    % omegaStep = 0.00001;
    %
    % omegaIndex = 1:length(w);
    omega = 0.00001:0.00001:1.6;

    k3p = sqrt(kp .^ 2 .* (NNp .^ 2 ./ omega .^ 2 - 1));
    mp = k3p * 1j;

    a11_mpIsRe = a11;
    a11_mpIsIm = 1 - a11_mpIsRe;
    a11 = a11_mpIsRe .* cosh(mp .* hp) + a11_mpIsIm .* cos(k3p .* hp);

    a12_mpIsRe = a12;
    a12_mpIsIm = 1 - a12_mpIsRe;
    a12 = a12_mpIsRe .* sinh(mp .* hp) ./ mp + a12_mpIsIm .* sin(k3p .* hp) ./ k3p;

    a21_mpIsRe = a21;
    a21_mpIsIm = 1 - a21_mpIsRe;
    a21 = a21_mpIsRe .* mp .* sinh(mp .* hp) - a21_mpIsIm .* sin(k3p .* hp) .* k3p;

    a22_mpIsRe = a22;
    a22_mpIsIm = 1 - a22_mpIsRe;
    a22 = a22_mpIsRe .* cosh(mp .* hp) + a22_mpIsIm .* cos(k3p .* hp);

    % a11(row_w_eq_Np, col_w_eq_Np) = 0;
    % a12(row_w_eq_Np, col_w_eq_Np) = 0;
    % a21(row_w_eq_Np, col_w_eq_Np) = 0;
    % a22(row_w_eq_Np, col_w_eq_Np) = 0;

    clear a11_mpIsRe a11_mpIsIm a12_mpIsRe a12_mpIsIm a21_mpIsRe a21_mpIsIm a22_mpIsRe a22_mpIsIm apJudge

    result11 = a11(1, :);
    result12 = a12(1, :);
    result21 = a21(1, :);
    result22 = a22(1, :);

    for i = 2:size(a11, 1)
        result11 = result11 .* a11(i, :) + result12 .* a21(i, :);
        result12 = result11 .* a12(i, :) + result12 .* a22(i, :);
        result21 = result21 .* a11(i, :) + result22 .* a21(i, :);
        result22 = result21 .* a12(i, :) + result22 .* a22(i, :);
    end

    f12(kp / 5, :) = result12;

end

toc
load F.mat
r1 = a11 - a11_test;
r2 = a12 - a12_test;
r3 = a21 - a21_test;
r4 = a22 - a22_test;
% function result12 = mulMatrix(a11, a12, a21, a22)
%     format long;
%     result11 = a11(1, :);
%     result12 = a12(1, :);
%     result21 = a21(1, :);
%     result22 = a22(1, :);

%     for i = 2:size(a11, 1)
%         result11 = result11 .* a11(i, :) + result12 .* a21(i, :);
%         result12 = result11 .* a12(i, :) + result12 .* a22(i, :);
%         result21 = result21 .* a11(i, :) + result22 .* a21(i, :);
%         result22 = result21 .* a12(i, :) + result22 .* a22(i, :);
%     end

% end
