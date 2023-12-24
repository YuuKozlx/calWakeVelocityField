clc;
clear;
close;
% load '../求解特征函数/求解所有的特征函数/eigenfunction.mat'
load eigenfunction_fit.mat

a = 0.20;
b = 0.05;
U = 0.90;
x_index = 1;
tic
for x_bar = 0:0.1:10
    
% mode 1
wm = eigenfunction(1).w;
cpm = eigenfunction(1).cp;
cgm = eigenfunction(1).cg;
Q = 4 * pi * b^2 * a * (sin(wm*a/U) - wm * a / U .* cos(wm*a/U)) ./ (wm * a / U).^3;

k = eigenfunction(1).k;

ky = sqrt(k.^2-(wm / U).^2);
kyisnotcomplex_index = (find(conj(ky)+ky))';
kyiscomplex_index = 1:kyisnotcomplex_index(1)-1;
k = k(kyisnotcomplex_index);
wm = wm(kyisnotcomplex_index);
cpm = cpm(kyisnotcomplex_index);
cgm = cgm(kyisnotcomplex_index);
Q = Q(kyisnotcomplex_index);
ky = sqrt(k.^2-(wm / U).^2);
ky_fit = 0:1:100;
vx_mode1 = zeros(801,length(kyisnotcomplex_index));
vy_mode1 = zeros(801,length(kyisnotcomplex_index));
vz_mode1 = zeros(801,length(kyisnotcomplex_index));
for i = 1:1:801
    vx_mode1(i,:) = 1i / (2 .* U^2) .* (Q .* cpm.^5 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(1).phi(i,kyisnotcomplex_index) .* eigenfunction(1).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
    vy_mode1(i,:) = 1i / (2 .* U).*sqrt(1./cpm.^2 - 1./U^2).* (Q .* cpm.^5 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(1).phi(i,kyisnotcomplex_index) .* eigenfunction(1).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
    vz_mode1(i,:) = -1 / (2 .* U) .* (Q .* cpm.^4 .* k.^2) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(1).phi(i,kyisnotcomplex_index) .* eigenfunction(1).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
end
% vx_mode1()
for i = 1:100:801
[fitobject_Re,~,~] = fit(ky,real(vx_mode1(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vx_mode1(i,:)'),'spline');
vx_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vx_mode_Im(i,:) = feval(fitobject_Im,ky_fit');

[fitobject_Re,~,~] = fit(ky,real(vy_mode1(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vy_mode1(i,:)'),'spline');
vy_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vy_mode_Im(i,:) = feval(fitobject_Im,ky_fit');

[fitobject_Re,~,~] = fit(ky,real(vz_mode1(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vz_mode1(i,:)'),'spline');
vz_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vz_mode_Im(i,:) = feval(fitobject_Im,ky_fit');
end
vx_mode1_fit = vx_mode_Re + 1i*vx_mode_Im;
vy_mode1_fit = vy_mode_Re + 1i*vy_mode_Im;
vz_mode1_fit = vz_mode_Re + 1i*vz_mode_Im;

%% mode 2---------------------------------------------------------------%
% 
wm = eigenfunction(2).w;
cpm = eigenfunction(2).cp;
cgm = eigenfunction(2).cg;
Q = 4 * pi * b^2 * a * (sin(wm*a/U) - wm * a / U .* cos(wm*a/U)) ./ (wm * a / U).^3;

k = eigenfunction(2).k;

ky = sqrt(k.^2-(wm / U).^2);
kyisnotcomplex_index = (find(conj(ky)+ky))';
kyiscomplex_index = 1:kyisnotcomplex_index(1)-1;
k = k(kyisnotcomplex_index);
wm = wm(kyisnotcomplex_index);
cpm = cpm(kyisnotcomplex_index);
cgm = cgm(kyisnotcomplex_index);
Q = Q(kyisnotcomplex_index);
ky = sqrt(k.^2-(wm / U).^2);
ky_fit = 0:1:100;
vx_mode2 = zeros(801,length(kyisnotcomplex_index));
vy_mode2 = zeros(801,length(kyisnotcomplex_index));
vz_mode2 = zeros(801,length(kyisnotcomplex_index));
for i = 1:1:801
    vx_mode2(i,:) = 1i / (2 .* U^2) .* (Q .* cpm.^5 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(2).phi(i,kyisnotcomplex_index) .* eigenfunction(2).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
    vy_mode2(i,:) = 1i / (2 .* U).*sqrt(1./cpm.^2 - 1./U^2).* (Q .* cpm.^5 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(2).phi(i,kyisnotcomplex_index) .* eigenfunction(2).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
    vz_mode2(i,:) = -1 / (2 .* U) .* (Q .* cpm.^4 .* k.^2) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(2).phi(i,kyisnotcomplex_index) .* eigenfunction(2).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
end
% vx_mode1()
for i = 1:100:801
[fitobject_Re,~,~] = fit(ky,real(vx_mode2(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vx_mode2(i,:)'),'spline');
vx_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vx_mode_Im(i,:) = feval(fitobject_Im,ky_fit');

[fitobject_Re,~,~] = fit(ky,real(vy_mode2(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vy_mode2(i,:)'),'spline');
vy_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vy_mode_Im(i,:) = feval(fitobject_Im,ky_fit');

[fitobject_Re,~,~] = fit(ky,real(vz_mode2(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vz_mode2(i,:)'),'spline');
vz_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vz_mode_Im(i,:) = feval(fitobject_Im,ky_fit');
end
vx_mode2_fit = vx_mode_Re + 1i*vx_mode_Im;
vy_mode2_fit = vy_mode_Re + 1i*vy_mode_Im;
vz_mode2_fit = vz_mode_Re + 1i*vz_mode_Im;
%%---------------------------------------------------------------%

%% mode 3---------------------------------------------------------------%
% 
wm = eigenfunction(3).w;
cpm = eigenfunction(3).cp;
cgm = eigenfunction(3).cg;
Q = 4 * pi * b^2 * a * (sin(wm*a/U) - wm * a / U .* cos(wm*a/U)) ./ (wm * a / U).^3;

k = eigenfunction(3).k;

ky = sqrt(k.^2-(wm / U).^2);
kyisnotcomplex_index = (find(conj(ky)+ky))';
kyiscomplex_index = 1:kyisnotcomplex_index(1)-1;
k = k(kyisnotcomplex_index);
wm = wm(kyisnotcomplex_index);
cpm = cpm(kyisnotcomplex_index);
cgm = cgm(kyisnotcomplex_index);
Q = Q(kyisnotcomplex_index);
ky = sqrt(k.^2-(wm / U).^2);
ky_fit = 0:1:100;
vx_mode3 = zeros(801,length(kyisnotcomplex_index));
vy_mode3 = zeros(801,length(kyisnotcomplex_index));
vz_mode3 = zeros(801,length(kyisnotcomplex_index));
for i = 1:1:801
    vx_mode3(i,:) = 1i / (2 .* U^2) .* (Q .* cpm.^5 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(3).phi(i,kyisnotcomplex_index) .* eigenfunction(3).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
    vy_mode3(i,:) = 1i / (2 .* U).*sqrt(1./cpm.^2 - 1./U^2).* (Q .* cpm.^5 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(3).phi(i,kyisnotcomplex_index) .* eigenfunction(3).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
    vz_mode3(i,:) = -1 / (2 .* U) .* (Q .* cpm.^4 .* k.^2) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(3).phi(i,kyisnotcomplex_index) .* eigenfunction(3).dphi(401,kyisnotcomplex_index))' .* exp(-1i.*wm.*x_bar/U);
end
% vx_mode1()
for i = 1:100:801
[fitobject_Re,~,~] = fit(ky,real(vx_mode3(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vx_mode3(i,:)'),'spline');
vx_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vx_mode_Im(i,:) = feval(fitobject_Im,ky_fit');

[fitobject_Re,~,~] = fit(ky,real(vy_mode3(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vy_mode3(i,:)'),'spline');
vy_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vy_mode_Im(i,:) = feval(fitobject_Im,ky_fit');

[fitobject_Re,~,~] = fit(ky,real(vz_mode3(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(vz_mode3(i,:)'),'spline');
vz_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
vz_mode_Im(i,:) = feval(fitobject_Im,ky_fit');
end
vx_mode3_fit = vx_mode_Re + 1i*vx_mode_Im;
vy_mode3_fit = vy_mode_Re + 1i*vy_mode_Im;
vz_mode3_fit = vz_mode_Re + 1i*vz_mode_Im;
%%---------------------------------------------------------------%
% mode 1
vx_test = vx_mode1_fit(1:100:801,:);
vx_amp = [conj(vx_test.');fliplr(vx_test.')]*0.5;
vx1 = ifft(vx_amp,'symmetric');

vy_test = vx_mode1_fit(1:100:801,:);
vy_amp = [conj(vy_test.');fliplr(vy_test.')]*0.5;
vy1 = ifft(vy_amp,'symmetric');

vz_test = vz_mode1_fit(1:100:801,:);
vz_amp = [conj(vz_test.');fliplr(vz_test.')]*0.5;
vz1 = ifft(vz_amp,'symmetric');
% mode 2
vx_test = vx_mode2_fit(1:100:801,:);
vx_amp = [conj(vx_test.');fliplr(vx_test.')]*0.5;
vx2 = ifft(vx_amp,'symmetric');

vy_test = vx_mode2_fit(1:100:801,:);
vy_amp = [conj(vy_test.');fliplr(vy_test.')]*0.5;
vy2 = ifft(vy_amp,'symmetric');

vz_test = vz_mode2_fit(1:100:801,:);
vz_amp = [conj(vz_test.');fliplr(vz_test.')]*0.5;
vz2 = ifft(vz_amp,'symmetric');

% mode 3
vx_test = vx_mode3_fit(1:100:801,:);
vx_amp = [conj(vx_test.');fliplr(vx_test.')]*0.5;
vx3 = ifft(vx_amp,'symmetric');

vy_test = vx_mode3_fit(1:100:801,:);
vy_amp = [conj(vy_test.');fliplr(vy_test.')]*0.5;
vy3 = ifft(vy_amp,'symmetric');

vz_test = vz_mode3_fit(1:100:801,:);
vz_amp = [conj(vz_test.');fliplr(vz_test.')]*0.5;
vz3 = ifft(vz_amp,'symmetric');

result_vx(x_index,:,:) = vx1 + vx2 + vx3;
result_vy(x_index,:,:) = vy1 + vy2 + vy3;
result_vz(x_index,:,:) = vz1 + vz2 + vz3;

x_index = x_index+1;

end
toc
x = 0:0.1:10;
result_all_vx = [fliplr(result_vx), result_vx];
result_all_vy = [-fliplr(result_vy), result_vy];
result_all_vz = [fliplr(result_vz), result_vz];
figure(1)
tiledlayout(2,2)
nexttile
surf(result_all_vx(:,100:300,2));shading interp;view([0 90]);
nexttile
surf(result_all_vy(:,100:300,2));shading interp;view([0 90]);
nexttile
surf(result_all_vz(:,100:300,2));shading interp;view([0 90]);