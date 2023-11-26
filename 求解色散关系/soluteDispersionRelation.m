clc; clear;
format long;
temp_load = load('../sigmoid函数产生NV频率/hp.mat'); % 量化后的NV频率数据对应的深度，对应于论文中的hp
hp = temp_load.hp;
temp_load = load('../sigmoid函数产生NV频率/Np.mat'); % 量化后的NV频率数据对应的Np，对应于论文中的Np
Np = temp_load.Np;

step_kp = 2.5;
step_omega = 0.00001;


% 采用并行计算
% 创建进程池
% poolobj = gcp('nocreate');
% 
% if isempty(poolobj)
%     parpool('local',4);
% end
tic
for kpIndex = 1:40

    for omegaIndex = 1:160000
        f12(kpIndex, omegaIndex) = cal_f12(kpIndex * step_kp, omegaIndex * step_omega, hp, Np);
    end

end

toc

function f12 = cal_f12(kp, omega, hp, Np)
    k3p = sqrt(kp .^ 2 .* (Np .^ 2 / omega .^ 2 - 1));
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

    f12 = A(1, 2);
end
