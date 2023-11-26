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

apJudge = sign(w_repmat .^ 2 - NNp * w);
apJudge(apJudge == 0) = 1;
%  w>Np, apJudge = 1 mp为实数
%  w<Np, apJudge = -1 mp为虚数
a11 = apJudge;
a12 = apJudge;
a21 = apJudge;
a22 = apJudge;

kp = 5;

% omegaStep = 0.00001;
%
% omegaIndex = 1:length(w);
omega = 0.00001:0.00001:1.6;

k3p = sqrt(kp .^ 2 .* (NNp .^ 2 ./ omega .^ 2 - 1));
mp = k3p * 1j;



a11_mpIsRe = a11;
a11_mpIsRe(a11_mpIsRe == -1) = 0;
a11_mpIsIm = a11 - a11_mpIsRe;
a11 = a11_mpIsRe .* cosh(mp .* hp) + a11_mpIsIm .* mp .* sinh(mp .* hp);

a12_mpIsRe = a12;
a12_mpIsRe(a12_mpIsRe == -1) = 0;
a12_mpIsIm = a12 - a12_mpIsRe;
a12 = a12_mpIsRe .* sinh(mp .* hp) ./ mp + a12_mpIsIm .* cosh(mp .* hp);

a21_mpIsRe = a21;
a21_mpIsRe(a21_mpIsRe == -1) = 0;
a21_mpIsIm = a21 - a21_mpIsRe;
a21 = a21_mpIsRe .* mp .* sinh(mp .* hp) + a21_mpIsIm .* cosh(mp .* hp);

a22_mpIsRe = a22;
a22_mpIsRe(a22_mpIsRe == -1) = 0;
a22_mpIsIm = a22 - a22_mpIsRe;
a22 = a22_mpIsRe .* cosh(mp .* hp) + a22_mpIsIm .* mp .* sinh(mp .* hp);


A(:, 1, :) = a11;
A(:, 2, :) = a12;
A(:, 3, :) = a21;
A(:, 4, :) = a22;


% a11(a11 == 1) = cosh(mp .* hp);
% a12(a12 == 1) = sinh(mp .* hp) ./ mp;
% a21(a21 == 1) = mp .* sinh(mp .* hp);
% a22(a22 == 1) = cosh(mp .* hp);
% 
% a11(a11 == -1) = cos(k3p .* hp);
% a12(a12 == -1) = sin(k3p .* hp) ./ k3p;
% a21(a21 == -1) = -k3p .* sin(k3p .* hp);
% a22(a22 == -1) = cos(k3p .* hp);
