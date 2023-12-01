clc;
clear;
load '..\sigmoid函数产生NV频率\hp.mat'
load '..\sigmoid函数产生NV频率\Np.mat'
load cp_k.mat
Np = Np';
hp = hp';

cp0 = cp(1, 2);
k = 0;
Nv = Np(1:15);
alpha = sqrt(abs(Nv .^ 2 / cp0 .^ 2 - k ^ 2));

r = zeros(15, 1);

for i = 1:length(hp)

    if (Nv(i) .^ 2 / cp0 ^ 2 - k ^ 2) > 0
        r(i) = 1i * alpha(i);
    elseif (Nv(i) .^ 2 / cp0 ^ 2 - k ^ 2) < 0
        r(i) = alpha(i);
    end

end

C = 0.001 * 1i;

s = eye(2) * [C; -C];

ApBp = zeros(2, 15);
ApBp(:, 1) = s;

for i = 1:14
    s = 0.5 * [(1 + r(i + 1) / r(i)) * exp(r(i) * hp(i)), (1 - r(i + 1) / r(i)) * exp(-r(i) * hp(i)); ...
                   (1 - r(i + 1) / r(i)) * exp(r(i) * hp(i)), (1 + r(i + 1) / r(i)) * exp(-r(i) * hp(i))] * s;
    ApBp(:, i + 1) = s;

end

phi_h15 = s(1, 1) * exp(r(15) * hp(15)) + s(2, 1) * exp(-r(15) * hp(15));

ApBp = transpose(ApBp);

z = (0:0.001:0.8)';
y = phi(z, ApBp, r);

% 0.001000 + 0.000000i	0.001000 + 0.000000i
% 0.001000 + 0.000301i	0.001000 - 0.000301i
% 0.000972 + 0.000575i	0.000972 - 0.000575i
% 0.000944 + 0.000825i	0.000944 - 0.000825i
% 0.000900 + 0.001091i	0.000900 - 0.001091i
% 0.000839 + 0.001366i	0.000839 - 0.001366i
% 0.000748 + 0.001655i	0.000748 - 0.001655i
% 0.000619 + 0.001951i    0.000619 - 0.001951i
% -0.000014 + 0.001791i	-0.000014 - 0.0017910i
% -0.000164 + 0.001529i	-0.000164 - 0.001529i
% -0.000264 + 0.001262i	-0.000264 - 0.001262i
% -0.000333 + 0.000997i	-0.000333 - 0.000997i
% -0.000380 + 0.000735i	-0.000380 - 0.000735i
% -0.000416 + 0.000477i	-0.000416 - 0.000477i
% -0.000453 + 0.000011i	-0.000453 - 0.000011i

function y = phi(z, ApBp, r)
    Ap = ApBp(:, 1);
    Bp = ApBp(:, 2);
    y = zeros(size(z));

    condition1 = (z >= 0 & z < 0.138);
    y(condition1) = Ap(1) * exp(r(1) * z(condition1)) + Bp(1) * exp(-r(1) * z(condition1));

    condition2 = (z >= 0.138 & z < 0.176);
    y(condition2) = Ap(2) * exp(r(2) * (z(condition2) - 0.138)) + Bp(2) * exp(-r(2) * (z(condition2) - 0.138));

    condition3 = (z >= 0.176 & z < 0.190);
    y(condition3) = Ap(3) * exp(r(3) * (z(condition3) - 0.176)) + Bp(3) * exp(-r(3) * (z(condition3) - 0.176));

    condition4 = (z >= 0.190 & z < 0.202);
    y(condition4) = Ap(4) * exp(r(4) * (z(condition4) - 0.190)) + Bp(4) * exp(-r(4) * (z(condition4) - 0.190));

    condition5 = (z >= 0.202 & z < 0.212);
    y(condition5) = Ap(5) * exp(r(5) * (z(condition5) - 0.202)) + Bp(5) * exp(-r(5) * (z(condition5) - 0.202));

    condition6 = (z >= 0.212 & z < 0.222);
    y(condition6) = Ap(6) * exp(r(6) * (z(condition6) - 0.212)) + Bp(6) * exp(-r(6) * (z(condition6) - 0.212));

    condition7 = (z >= 0.222 & z < 0.232);
    y(condition7) = Ap(7) * exp(r(7) * (z(condition7) - 0.222)) + Bp(7) * exp(-r(7) * (z(condition7) - 0.222));

    condition8 = (z >= 0.232 & z < 0.268);
    y(condition8) = Ap(8) * exp(r(8) * (z(condition8) - 0.232)) + Bp(8) * exp(-r(8) * (z(condition8) - 0.232));

    condition9 = (z >= 0.268 & z < 0.279);
    y(condition9) = Ap(9) * exp(r(9) * (z(condition9) - 0.268)) + Bp(9) * exp(-r(9) * (z(condition9) - 0.268));

    condition10 = (z >= 0.279 & z < 0.289);
    y(condition10) = Ap(10) * exp(r(10) * (z(condition10) - 0.279)) + Bp(10) * exp(-r(10) * (z(condition10) - 0.279));

    condition11 = (z >= 0.289 & z < 0.299);
    y(condition11) = Ap(11) * exp(r(11) * (z(condition11) - 0.289)) + Bp(11) * exp(-r(11) * (z(condition11) - 0.289));

    condition12 = (z >= 0.299 & z < 0.310);
    y(condition12) = Ap(12) * exp(r(12) * (z(condition12) - 0.299)) + Bp(12) * exp(-r(12) * (z(condition12) - 0.299));

    condition13 = (z >= 0.310 & z < 0.325);
    y(condition13) = Ap(13) * exp(r(13) * (z(condition13) - 0.310)) + Bp(13) * exp(-r(13) * (z(condition13) - 0.310));

    condition14 = (z >= 0.325 & z < 0.362);
    y(condition14) = Ap(14) * exp(r(14) * (z(condition14) - 0.325)) + Bp(14) * exp(-r(14) * (z(condition14) - 0.325));

    condition15 = (z >= 0.362 & z < 0.800);
    y(condition15) = Ap(15) * exp(r(15) * (z(condition15) - 0.362)) + Bp(15) * exp(-r(15) * (z(condition15) - 0.362));

end
