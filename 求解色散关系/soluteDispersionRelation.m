clc; clear;
format long;
temp_load = load('../sigmoid函数产生NV频率/hp.mat'); % 量化后的NV频率数据对应的深度，对应于论文中的hp
hp = temp_load.hp;
temp_load = load('../sigmoid函数产生NV频率/Np.mat'); % 量化后的NV频率数据对应的Np，对应于论文中的Np
Np = temp_load.Np;

syms kp;
syms omega;
%  NV频率

kpIndex = 1;
omegaIndex = 1;

for kp = 0.001:0.001:100

    for omega = 0.01:0.01:1.8
        k3p = kp .^ 2 .* abs(Np .^ 2 / omega .^ 2 - 1);
        mp = k3p * 1j;

        for index = 1:length(hp)

            if omega > Np(index)
                a11(index, 1) = cosh(mp(1, index) * hp(index));
                a12(index, 1) = sinh(mp(1, index) * hp(index)) ./ mp(1, index);
                a21(index, 1) = mp(index) * sinh(mp(1, index) * hp(index));
                a22(index, 1) = cosh(mp(1, index) * hp(index));
            elseif omega < Np(index)
                a11(index, 1) = cos(k3p(1, index) * hp(index));
                a12(index, 1) = sin(k3p(1, index) * hp(index)) / k3p(1, index);
                a21(index, 1) = -k3p(1, index) * sin(k3p(1, index) * hp(index));
                a22(index, 1) = cos(k3p(1, index) * hp(index));
            end

        end

        A = [a11(1, 1), a12(1, 1); a21(1, 1), a22(1, 1)];

        for index = 2:length(hp)
            A = [a11(index, 1), a12(index, 1); a21(index, 1), a22(index, 1)] * A;
        end

        f12(kpIndex,omegaIndex) = abs(A(1, 2));
        omegaIndex = omegaIndex + 1;

    end
    kpIndex = kpIndex + 1;
    omegaIndex = 1;
end
