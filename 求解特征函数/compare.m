clc;clear;
q_0 = 0.001;
q_end = 0.001;
k_index = 1;
mode = 3;
run test3_TH.m
run test4_TH.m
load phi_0.mat
load phi_end.mat

q_end = sign(phi_0(400)*phi_end(400))*q_end*max(abs(phi_0(300:500)))/max(abs(phi_end(300:500)));
run test3_TH.m
run test4_TH.m
load phi_0.mat
load phi_end.mat

phi = [phi_0,phi_end];
figure(5)
plot(phi);
xlim([0 800])
legend('phi 0','phi end','Location','northeastoutside');

figure(6)
phi_std(1:350) = phi_0(1:350);
phi_std(351:451) = (phi_0(351:451)+phi_end(351:451))/2;
phi_std(452:801) = phi_end(452:end);
plot(phi_std);
xlim([0 800])