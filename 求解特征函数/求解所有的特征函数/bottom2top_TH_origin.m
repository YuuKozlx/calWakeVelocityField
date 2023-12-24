% TH法求解，从底部开始算到顶部

load '..\..\函数产生NV频率\hp.mat'
load '..\..\函数产生NV频率\Np.mat'
load '..\..\求解色散关系\色散关系 ω_k 结果\cp.mat'
load '..\..\求解色散关系\色散关系 ω_k 结果\cg.mat'

Np = Np';
hp = hp';


cp0 = cp(k_index, mode+1);
k = k_index * 1 - 0.95;
Np = Np(1:end-1);

alpha = sqrt(abs(Np.^2/cp0.^2-k^2));

zp = zeros(length(Np), 1);

for i = 1:length(hp)
    zp(i, :) = sum(hp(1:i));
end

mp = zeros(length(Np), 1);
rp = zeros(length(Np), 1);

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
q = q_end;
% q = 0.001;

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




for i = length(Np) - 1:-1:1

    if imag(mp(i)) == 0
        yp0(:,i) = [y{i + 1, 1}(q, 0); dy{i + 1, 1}(q, 0)];
        temp = 0.5 * [exp(-mp(i).*hp(i)), 1 / mp(i) .* exp(-mp(i).*hp(i)); exp(mp(i).*hp(i)), -1 / mp(i) .* exp(mp(i).*hp(i))] * [y{i + 1, 1}(q, 0); dy{i + 1, 1}(q, 0)];
        Ap(i, :) = temp(1);
        Bp(i, :) = temp(2);
        y{i, 1} = @(An, z) Ap(i) .* exp(mp(i).*z) + Bp(i) .* exp(-mp(i).*z);
        dy{i, 1} = @(An, z) mp(i) .* Ap(i) .* exp(mp(i).*z) - mp(i) .* Bp(i) .* exp(-mp(i).*z);
    elseif imag(mp(i)) ~= 0
        yp0(:,i) = [y{i + 1, 1}(q, 0); dy{i + 1, 1}(q, 0)];
        temp = [sin(rp(i).*hp(i)), 1 / rp(i) .* cos(rp(i).*hp(i)); cos(rp(i).*hp(i)), -1 / rp(i) .* sin(rp(i).*hp(i))] * [y{i + 1, 1}(q, 0); dy{i + 1, 1}(q, 0)];
        Ap(i, :) = temp(1);
        Bp(i, :) = temp(2);

        y{i, 1} = @(An, z) Ap(i) .* sin(rp(i).*z) + Bp(i) .* cos(rp(i).*z);
        dy{i, 1} = @(An, z) rp(i) .* Ap(i) .* cos(rp(i).*z) - rp(i) .* Bp(i) .* sin(rp(i).*z);
    end

end



phi = [];
dphi = [];
z = (0:0.001:0.8)';
for i = 1:length(Np)
    if i == 1
        temp = y{i,1}(q,z(z >= 0 & z < zp(1)));
        dtemp = dy{i,1}(q,z(z >= 0 & z < zp(1)));
    elseif i == length(Np)
        temp = y{i,1}(q,z(z >= zp(i-1) & z <= zp(i)) - zp(i-1));
        dtemp = dy{i,1}(q,z(z >= zp(i-1) & z <= zp(i)) - zp(i-1));
    else
        temp = y{i,1}(q,z(z >= zp(i-1) & z < zp(i)) - zp(i-1));
        dtemp = dy{i,1}(q,z(z >= zp(i-1) & z < zp(i)) - zp(i-1));
    end
    phi = [phi;temp];
    dphi = [dphi;dtemp];
end



phi_end = phi;
dphi_end = dphi;
save('phi_end.mat',"phi_end");
save('dphi_end.mat',"dphi_end");