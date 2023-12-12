clc;
clear;
load '..\sigmoid函数产生NV频率\hp.mat'
load '..\sigmoid函数产生NV频率\Np.mat'
load '..\求解色散关系\色散关系 ω_k 结果\cp.mat'
Np = Np';
hp = hp';

k_index = 11;
cp0 = cp(k_index, 2);
k = (k_index - 1) * 0.5;
Np = Np(1:15);

alpha = sqrt(abs(Np.^2/cp0.^2-k^2));

zp = zeros(15, 1);

for i = 1:length(hp)
    zp(i, :) = sum(hp(1:i));
end

mp = zeros(15, 1);
rp = zeros(15, 1);

for i = 1:length(hp)

    if (Np(i).^2 / cp0^2 - k^2) > 0
        mp(i) = 1i * alpha(i);
        rp(i) = -1i * mp(i);
    elseif (Np(i).^2 / cp0^2 - k^2) < 0
        mp(i) = alpha(i);
        rp(i) = mp(i);
    end

end

syms An;
syms z;
q = -1;

if imag(mp(end)) == 0
    Bn = @(An) -An * exp(2*mp(end)*hp(end));
elseif imag(mp(end)) ~= 0
    Bn = @(An) -An * tan(rp(end)*hp(end));
end

if imag(mp(end)) == 0
    y{length(Np), 1} = @(An, z) An .* exp(mp(end).*z) + Bn(An) .* exp(-mp(end).*z);
    dy{length(Np), 1} = @(An, z) mp(end) .* An .* exp(mp(end).*z) - mp(end) .* Bn(An) .* exp(-mp(end).*z);
elseif imag(mp(end)) ~= 0
    y{length(Np), 1} = @(An, z) An .* sin(rp(end).*z) + Bn(An) .* cos(rp(end).*z);
    dy{length(Np), 1} = @(An, z) rp(end) .* An .* cos(rp(end).*z) - rp(end) .* Bn(An) .* sin(rp(end).*z);
end

i = 1;

for z = 0.001:0.001:0.438
    phi(i) = y{15, 1}(1, z);
    dphi(i) = dy{15, 1}(1, z);
    i = i + 1;
end

for i = length(Np) - 1:-1:1

    if imag(mp(i)) == 0
        temp = 0.5 * [exp(-mp(i).*hp(i)), 1 / mp(i) .* exp(-mp(i).*hp(i)); mp(i) .* exp(mp(i).*hp(i)), -1 / mp(i) .* exp(-mp(i).*hp(i))] * [y{i + 1, 1}(q, 0); dy{i + 1, 1}(q, 0)];
        Ap(i, :) = temp(1);
        Bp(i, :) = temp(2);
        y{i, 1} = @(An, z) Ap(i) .* exp(mp(i).*z) + Bp(i) .* exp(-mp(i).*z);
        dy{i, 1} = @(An, z) mp(i) .* Ap(i) .* exp(mp(i).*z) - mp(i) .* Bp(i) .* exp(-mp(i).*z);
    elseif imag(mp(i)) ~= 0
        temp = [sin(rp(i).*hp(i)), 1 / rp(i) .* cos(rp(i).*hp(i)); cos(rp(i).*hp(i)), -1 / rp(i) .* sin(rp(i).*hp(i))] * [y{i + 1, 1}(q, 0); dy{i + 1, 1}(q, 0)];
        Ap(i, :) = temp(1);
        Bp(i, :) = temp(2);

        y{i, 1} = @(An, z) Ap(i) .* sin(rp(i).*z) + Bp(i) .* cos(rp(i).*z);
        dy{i, 1} = @(An, z) rp(i) .* Ap(i) .* cos(rp(i).*z) - rp(i) .* Bp(i) .* sin(rp(i).*z);
    end

end
phi = [];
dphi = [];
for j = 1:length(hp)
    if j == 1
        z = 0:0.0001:hp(j);
    else
        z = 0.0001:0.0001:hp(j);
    end
    z = z';
    phi = [phi; y{i, 1}(1, z)];
    dphi = [dphi; dy{i, 1}(1, z)];
end

a1 = y{14, 1}(1, 0);
a2 = dy{14, 1}(1, 0);
b1 = y{13, 1}(1, hp(13));
b2 = dy{13, 1}(1, hp(13));

a1 - b1
a2 - b2
