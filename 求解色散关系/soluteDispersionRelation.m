clc; clear;
format long;
temp_load = load('../sigmoid函数产生NV频率/hp.mat'); % 量化后的NV频率数据对应的深度，对应于论文中的hp
hp = temp_load.hp;
temp_load = load('../sigmoid函数产生NV频率/Np.mat'); % 量化后的NV频率数据对应的Np，对应于论文中的Np
Np = temp_load.Np;
NNp = Np(1:15); % 只取前15个Np值。

f12 = zeros(20, 160000);
a11 = zeros(length(hp), 1);
a12 = zeros(length(hp), 1);
a21 = zeros(length(hp), 1);
a22 = zeros(length(hp), 1);

kpIndex = 1;
omegaIndex = 1;
tic

for kp = 5:5:100
    % kp = 80;
    for omega = 0.000005:0.00001:1.6
        k3p = sqrt(kp .^ 2 .* (NNp .^ 2 / omega .^ 2 - 1));
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

            % 注意：这种针对A的计算方法，当omega == Np时，理应有a11 = a12 = a21 = a22 = 0,如此在计算
            % F = a0*a1*a2...an-1*an时会有一个矩阵ai等于0矩阵，但是原论文中并没有对这一点进行处理
            % 为避免在计算中出现这种问题，可考虑在并行优化时给omega增加一个小于步长的偏置量，并在避免nv频率的取值与omega的各个取值重合。
            % 而程序中在遇到 omega == Np位置时作为替代，a11 a12 a21 a22 取值为前一个Omega 对应取值
            % 例如当 Np = 1，也就是 index =5 时 omega 若采用 0.00001：0.00001：1.6，则当omega取到1.00000时
            % 理应有a11(5，1) = a12(5，1)= a21(5，1)= a22(5，1) =0；但是判断条件时，跳过该步骤。因此这个对应取值为
            % omega = 0.99999时对应的 a11(5，1) a12(5，1) a21(5，1) a22(5，1).
        end

        A = [a11(1, 1), a12(1, 1); a21(1, 1), a22(1, 1)];

        for index = 2:length(hp)
            A = [a11(index, 1), a12(index, 1); a21(index, 1), a22(index, 1)] * A;
        end

        f12(kpIndex, omegaIndex) = A(1, 2);
        omegaIndex = omegaIndex + 1;

    end

    kpIndex = kpIndex + 1;
    omegaIndex = 1;
end

toc
