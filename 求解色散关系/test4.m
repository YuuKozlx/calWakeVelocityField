clc; clear;
format long;
temp_load = load('../sigmoid函数产生NV频率/hp.mat'); % 量化后的NV频率数据对应的深度，对应于论文中的hp
hp = temp_load.hp;
temp_load = load('../sigmoid函数产生NV频率/Np.mat'); % 量化后的NV频率数据对应的Np，对应于论文中的Np
Np = temp_load.Np;
NNp = Np(1:15);
a11 = zeros(15, 160000);
a12 = zeros(15, 160000);
a21 = zeros(15, 160000);
a22 = zeros(15, 160000);
f12 = zeros(20, 160000);

kpIndex = 1;
omegaIndex = 1;
tic
for kp = 5:5:100

for omega = 0.000005:0.00001:1.6
k3p = sqrt(kp .^ 2 .* (NNp .^ 2 / omega .^ 2 - 1));
mp = k3p * 1j;

for index = 1:length(hp)

    if omega > Np(index)
        a11(index, omegaIndex) = cosh(mp(1, index) * hp(index));
        a12(index, omegaIndex) = sinh(mp(1, index) * hp(index)) ./ mp(1, index);
        a21(index, omegaIndex) = mp(index) * sinh(mp(1, index) * hp(index));
        a22(index, omegaIndex) = cosh(mp(1, index) * hp(index));
    elseif omega < Np(index)
        a11(index, omegaIndex) = cos(k3p(1, index) * hp(index));
        a12(index, omegaIndex) = sin(k3p(1, index) * hp(index)) / k3p(1, index);
        a21(index, omegaIndex) = -k3p(1, index) * sin(k3p(1, index) * hp(index));
        a22(index, omegaIndex) = cos(k3p(1, index) * hp(index));
    end

end

A = [a11(1, omegaIndex), a12(1, omegaIndex); a21(1, omegaIndex), a22(1, omegaIndex)];

for index = 2:length(hp)
    A = [a11(index, omegaIndex), a12(index, omegaIndex); a21(index, omegaIndex), a22(index, omegaIndex)] * A;
end

f12(kpIndex, omegaIndex) = A(1, 2);
omegaIndex = omegaIndex + 1;

end
kpIndex = kpIndex + 1;
omegaIndex = 1;
end
toc
