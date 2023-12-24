clc;
clear;

depth_maxNp = 300;
depth_factor = 70;
for mode = 1:3
   
    for k_index = 1:100

        q_0 = 0.001;
        q_end = 0.001;
        run bottom2top_TH.m
        run top2bottom_TH.m
        load phi_0.mat
        load phi_end.mat
        load dphi_0.mat
        load dphi_end.mat

        q_end = sign(phi_0(depth_maxNp) * phi_end(depth_maxNp)) * q_end * max(abs(phi_0(depth_maxNp-100:depth_maxNp+100))) / max(abs(phi_end(depth_maxNp-100:depth_maxNp+100)));
        run bottom2top_TH.m
        run top2bottom_TH.m
        load phi_0.mat
        load phi_end.mat
        load dphi_0.mat
        load dphi_end.mat

       

        if mode == 1
            phi_std_mode_1(1:depth_maxNp - depth_factor) = phi_0(1:depth_maxNp - depth_factor);
            phi_std_mode_1(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) = (phi_0(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) + phi_end(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1)) / 2;
            phi_std_mode_1(depth_maxNp + depth_factor + 2:801) = phi_end(depth_maxNp + depth_factor + 2:end);
            phi_mode_1(:, k_index) = phi_std_mode_1';

            dphi_std_mode_1(1:depth_maxNp - depth_factor) = dphi_0(1:depth_maxNp - depth_factor);
            dphi_std_mode_1(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) = (dphi_0(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) + dphi_end(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1)) / 2;
            dphi_std_mode_1(depth_maxNp + depth_factor + 2:801) = dphi_end(depth_maxNp + depth_factor + 2:end);
            dphi_mode_1(:, k_index) = dphi_std_mode_1';
        elseif mode == 2
            phi_std_mode_2(1:depth_maxNp - depth_factor) = phi_0(1:depth_maxNp - depth_factor);
            phi_std_mode_2(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) = (phi_0(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) + phi_end(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1)) / 2;
            phi_std_mode_2(depth_maxNp + depth_factor + 2:801) = phi_end(depth_maxNp + depth_factor + 2:end);
            phi_mode_2(:, k_index) = phi_std_mode_2';

            dphi_std_mode_2(1:depth_maxNp - depth_factor) = dphi_0(1:depth_maxNp - depth_factor);
            dphi_std_mode_2(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) = (dphi_0(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) + dphi_end(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1)) / 2;
            dphi_std_mode_2(depth_maxNp + depth_factor + 2:801) = dphi_end(depth_maxNp + depth_factor + 2:end);
            dphi_mode_2(:, k_index) = dphi_std_mode_2';
        else
            phi_std_mode_3(1:depth_maxNp - depth_factor) = phi_0(1:depth_maxNp - depth_factor);
            phi_std_mode_3(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) = (phi_0(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) + phi_end(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1)) / 2;
            phi_std_mode_3(depth_maxNp + depth_factor + 2:801) = phi_end(depth_maxNp + depth_factor + 2:end);
            phi_mode_3(:, k_index) = phi_std_mode_3';

            dphi_std_mode_3(1:depth_maxNp - depth_factor) = dphi_0(1:depth_maxNp - depth_factor);
            dphi_std_mode_3(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) = (dphi_0(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1) + dphi_end(depth_maxNp - depth_factor + 1:depth_maxNp + depth_factor + 1)) / 2;
            dphi_std_mode_3(depth_maxNp + depth_factor + 2:801) = dphi_end(depth_maxNp + depth_factor + 2:end);
            dphi_mode_3(:, k_index) = dphi_std_mode_3';
        end

    end

end

% 求归一化后的特征函数
load '..\..\函数产生NV频率\nvQuant.mat';
load '..\..\求解色散关系\色散关系 ω_k 结果\w.mat';
nvQuant = nvQuant';

% mode 1
norm_coef = nvQuant .^ 2 .* phi_mode_1 .^ 2 * 0.001;
norm_coef = sum(norm_coef);
norm_coef = sqrt(norm_coef);
phi_mode_1 = phi_mode_1 ./ norm_coef;
dphi_mode_1 = dphi_mode_1 ./ norm_coef;

eigenfunction(1).mode = 'mode 1';
eigenfunction(1).k = cp(:, 1);
eigenfunction(1).w = w(1:end-1, 1);
eigenfunction(1).cp = cp(:, 2);
eigenfunction(1).cg = cg(:, 2);
eigenfunction(1).z = 0:0.001:0.800;
eigenfunction(1).phi = phi_mode_1;
eigenfunction(1).dphi = dphi_mode_1;

% mode 2
norm_coef = nvQuant .^ 2 .* phi_mode_2 .^ 2 * 0.001;
norm_coef = sum(norm_coef);
norm_coef = sqrt(norm_coef);
phi_mode_2 = phi_mode_2 ./ norm_coef;
dphi_mode_2 = dphi_mode_2 ./ norm_coef;

eigenfunction(2).mode = 'mode 2';
eigenfunction(2).k = cp(:, 1);
eigenfunction(2).w = w(1:end-1, 2);
eigenfunction(2).cp = cp(:, 3);
eigenfunction(2).cg = cg(:, 3);
eigenfunction(2).z = 0:0.001:0.800;
eigenfunction(2).phi = phi_mode_2;
eigenfunction(2).dphi = dphi_mode_2;

% mode 3

norm_coef = nvQuant .^ 2 .* phi_mode_3 .^ 2 * 0.001;
norm_coef = sum(norm_coef);
norm_coef = sqrt(norm_coef);
phi_mode_3 = phi_mode_3 ./ norm_coef;
dphi_mode_3 = dphi_mode_3 ./ norm_coef;

eigenfunction(3).mode = 'mode 3';
eigenfunction(3).k = cp(:, 1);
eigenfunction(3).w = w(1:end-1, 3);
eigenfunction(3).cp = cp(:, 4);
eigenfunction(3).cg = cg(:, 4);
eigenfunction(3).z = 0:0.001:0.800;
eigenfunction(3).phi = phi_mode_3;
eigenfunction(3).dphi = dphi_mode_3;

save('eigenfunction_origin.mat', 'eigenfunction');
