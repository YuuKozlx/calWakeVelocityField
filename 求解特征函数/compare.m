clc; clear;
q_0 = 0.001;
q_end = 0.001;
k_index = 198;
mode = 1;
run test3_TH.m
run test4_TH.m
load phi_0.mat
load phi_end.mat

q_end = sign(phi_0(300) * phi_end(300)) * q_end * max(abs(phi_0(200:400))) / max(abs(phi_end(200:400)));
run test3_TH.m
run test4_TH.m
load phi_0.mat
load phi_end.mat

phi = [phi_0, phi_end];
figure(5)
plot(phi);
xlim([0 800])
legend('phi 0', 'phi end', 'Location', 'northeastoutside');


depth_maxNp = 300;
depth_factor = 40;

figure(6)
phi_std(1:depth_maxNp-depth_factor) = phi_0(1:depth_maxNp-depth_factor);
phi_std(depth_maxNp-depth_factor+1:depth_maxNp+depth_factor+1) = (phi_0(depth_maxNp-depth_factor+1:depth_maxNp+depth_factor+1) + phi_end(depth_maxNp-depth_factor+1:depth_maxNp+depth_factor+1)) / 2;
phi_std(depth_maxNp-depth_factor+2:801) = phi_end(depth_maxNp-depth_factor+2:end);
plot(phi_std);
xlim([0 800])

figure(7)
dphi_std(1:depth_maxNp-depth_factor) = dphi_0(1:depth_maxNp-depth_factor);
dphi_std(depth_maxNp-depth_factor+1:depth_maxNp+depth_factor+1) = (dphi_0(depth_maxNp-depth_factor+1:depth_maxNp+depth_factor+1) + dphi_end(depth_maxNp-depth_factor+1:depth_maxNp+depth_factor+1)) / 2;
dphi_std(depth_maxNp-depth_factor+2:801) = dphi_end(depth_maxNp-depth_factor+2:end);
plot(dphi_std);
xlim([0 800])