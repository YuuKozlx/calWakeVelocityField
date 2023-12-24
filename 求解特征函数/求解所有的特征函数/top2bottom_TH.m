% TH法求解，从顶部开始算到底部
% clc;
% clear;
load '..\..\函数产生NV频率\hp.mat'
load '..\..\函数产生NV频率\Np.mat'
load '..\..\计算尾流速度场\cp.mat'
load '..\..\计算尾流速度场\cp.mat'
% load '..\..\求解色散关系\色散关系 ω_k 结果\cp.mat'
% load '..\..\求解色散关系\色散关系 ω_k 结果\cg.mat'
Np = Np';
hp = hp';


cp0 = cp(k_index, mode + 1);
% k = k_index * 1-0.95;
k = k_index * 0.05;
Np = Np(1:end - 1);

alpha = sqrt(abs(Np .^ 2 / cp0 .^ 2 - k ^ 2));

zp = zeros(length(Np), 1);

for i = 1:length(hp)
    zp(i, :) = sum(hp(1:i));
end

mp = zeros(length(Np), 1);
rp = zeros(length(Np), 1);

for i = 1:length(hp)

    if (Np(i) .^ 2 / cp0 ^ 2 - k ^ 2) > 0
        mp(i) = 1i * alpha(i);
        rp(i) = -1i * mp(i);
    elseif (Np(i) .^ 2 / cp0 ^ 2 - k ^ 2) < 0
        mp(i) = alpha(i);
        rp(i) = mp(i);
    end

end

syms A1;
syms z;
q = q_0;
% q = 0.001;

if imag(mp(1)) == 0
    B1 = @(A1) -A1;
elseif imag(mp(1)) ~= 0
    B1 = @(A1) 0;
end

if imag(mp(1)) == 0
    y{1, 1} = @(A1, z) A1 .* exp(mp(1) .* z) + B1(A1) .* exp(-mp(1) .* z);
    dy{1, 1} = @(A1, z) mp(1) .* A1 .* exp(mp(1) .* z) - mp(1) .* B1(A1) .* exp(-mp(1) .* z);
elseif imag(mp(1)) ~= 0
    y{1, 1} = @(A1, z) A1 .* sin(rp(1) .* z) + B1(A1) .* cos(rp(1) .* z);
    dy{1, 1} = @(A1, z) rp(1) .* A1 .* cos(rp(1) .* z) - rp(1) .* B1(A1) .* sin(rp(1) .* z);
end


for i = 2:1:length(Np)

    if imag(mp(i)) == 0

        yp_hp(:, i) = [y{i - 1, 1}(q, hp(i - 1)); dy{i - 1, 1}(q, hp(i - 1))];
        temp = 0.5 * [1, 1 / mp(i); 1, -1 / mp(i)] * [y{i - 1, 1}(q, hp(i - 1)); dy{i - 1, 1}(q, hp(i - 1))];
        Ap(i, :) = temp(1);
        Bp(i, :) = temp(2);
        y{i, 1} = @(A1, z) Ap(i) .* exp(mp(i) .* z) + Bp(i) .* exp(-mp(i) .* z);
        dy{i, 1} = @(A1, z) mp(i) .* Ap(i) .* exp(mp(i) .* z) - mp(i) .* Bp(i) .* exp(-mp(i) .* z);

    elseif imag(mp(i)) ~= 0

        yp_hp(:, i) = [y{i - 1, 1}(q, hp(i - 1)); dy{i - 1, 1}(q, hp(i - 1))];
        temp = [0, 1 / rp(i); 1, 0] * [y{i - 1, 1}(q, hp(i - 1)); dy{i - 1, 1}(q, hp(i - 1))];
        Ap(i, :) = temp(1);
        Bp(i, :) = temp(2);
        y{i, 1} = @(A1, z) Ap(i) .* sin(rp(i) .* z) + Bp(i) .* cos(rp(i) .* z);
        dy{i, 1} = @(A1, z) rp(i) .* Ap(i) .* cos(rp(i) .* z) - rp(i) .* Bp(i) .* sin(rp(i) .* z);

    end

end

% phi = [];
% dphi = [];
% for i = 1:length(hp)
%     if i == 1
%         z = (0:0.001:hp(i))';
%     else
%         z = (0.001:0.001:hp(i))';
%     end
%     temp = y{i,1}(q,z);
%     dtemp = dy{i,1}(q,z);
%     phi = [phi;temp];
%     dphi = [dphi;temp];
% end

phi = [];
dphi = [];
z = (0:0.001:0.8)';

for i = 1:length(Np)

    if i == 1
        temp = y{i, 1}(q, z(z >= 0 & z < zp(1)));
        dtemp = dy{i, 1}(q, z(z >= 0 & z < zp(1)));
    elseif i == length(Np)
        temp = y{i, 1}(q, z(z >= zp(i - 1) & z <= zp(i)) - zp(i - 1));
        dtemp = dy{i, 1}(q, z(z >= zp(i - 1) & z <= zp(i)) - zp(i - 1));
    else
        temp = y{i, 1}(q, z(z >= zp(i - 1) & z < zp(i)) - zp(i - 1));
        dtemp = dy{i, 1}(q, z(z >= zp(i - 1) & z < zp(i)) - zp(i - 1));
    end

    phi = [phi; temp];
    dphi = [dphi; dtemp];
end

phi_0 = phi;
dphi_0 = dphi;
save("phi_0.mat", "phi_0");
save("dphi_0.mat", "dphi_0");
