clc;
clear;
close;
% load '../求解特征函数/求解所有的特征函数/eigenfunction.mat'
load eigenfunction_fit.mat

a = 0.20;
b = 0.05;
U = 0.1;
x_index = 1;
for x_bar = 0:0.1:10
tic
% mode 1
wm = eigenfunction(1).w;
cpm = eigenfunction(1).cp;
cgm = eigenfunction(1).cg;
Q = 4 * pi * b^2 * a * (sin(wm*a/U) - wm * a / U .* cos(wm*a/U)) ./ (wm * a / U).^3;

k = (0.05:0.05:100)';

ky = sqrt(k.^2-(wm / U).^2);
ky_fit = 0.05:0.1:100;


eta_mode1 = zeros(801,2000);
for i = 1:100:801
    eta_mode1(i,:) = 1i / (2 .* U) .* (Q .* cpm.^3 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(1).phi(i,:) .* eigenfunction(1).dphi(401,:))' .* exp(-1i.*wm.*x_bar/U);
end
for i = 1:100:801
[fitobject_Re,~,~] = fit(ky,real(eta_mode1(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(eta_mode1(i,:)'),'spline');
eta_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
eta_mode_Im(i,:) = feval(fitobject_Im,ky_fit');
end

eta_mode1_fit = eta_mode_Re + 1i*eta_mode_Im;



% mode 2
wm = eigenfunction(2).w;
cpm = eigenfunction(2).cp;
cgm = eigenfunction(2).cg;
Q = 4 * pi * b^2 * a * (sin(wm*a/U) - wm * a / U .* cos(wm*a/U)) ./ (wm * a / U).^3;

k = (0.05:0.05:100)';

ky = sqrt(k.^2-(wm / U).^2);
ky_fit = 0.05:0.1:100;



eta_mode2 = zeros(801,2000);
for i = 1:100:801
    eta_mode2(i,:) = 1i / (2 .* U) .* (Q .* cpm.^3 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(2).phi(i,:) .* eigenfunction(2).dphi(401,:))' .* exp(-1i.*wm.*x_bar/U);
end
for i = 1:100:801
[fitobject_Re,~,~] = fit(ky,real(eta_mode2(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(eta_mode2(i,:)'),'spline');
eta_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
eta_mode_Im(i,:) = feval(fitobject_Im,ky_fit');
end

eta_mode2_fit = eta_mode_Re + 1i*eta_mode_Im;

% mode 3
wm = eigenfunction(3).w;
cpm = eigenfunction(3).cp;
cgm = eigenfunction(3).cg;
Q = 4 * pi * b^2 * a * (sin(wm*a/U) - wm * a / U .* cos(wm*a/U)) ./ (wm * a / U).^3;

k = (0.05:0.05:100)';

ky = sqrt(k.^2-(wm / U).^2);
ky_fit = 0.05:0.1:100;

eta_mode3 = zeros(801,2000);
for i = 1:100:801
    eta_mode3(i,:) = 1i / (2 .* U) .* (Q .* cpm.^3 .* k) ./ (1 - cpm .* cgm / U^2) .* (eigenfunction(3).phi(i,:) .* eigenfunction(3).dphi(401,:))' .* exp(-1i.*wm.*x_bar/U);
end
for i = 1:100:801
[fitobject_Re,~,~] = fit(ky,real(eta_mode3(i,:)'),'spline');
[fitobject_Im,~,~] = fit(ky,imag(eta_mode3(i,:)'),'spline');
eta_mode_Re(i,:) = feval(fitobject_Re,ky_fit');
eta_mode_Im(i,:) = feval(fitobject_Im,ky_fit');
end

eta_mode3_fit = eta_mode_Re + 1i*eta_mode_Im;

eta_test1 = eta_mode1_fit(401,:);
eta_amp = [conj(eta_test1),fliplr(eta_test1)]*0.5;
eta1 = ifft(eta_amp,'symmetric');

eta_test2 = eta_mode2_fit(401,:);
eta_amp = [conj(eta_test1),fliplr(eta_test1)]*0.5;
eta2 = ifft(eta_amp,'symmetric');

eta_test3 = eta_mode3_fit(401,:);
eta_amp = [conj(eta_test1),fliplr(eta_test1)]*0.5;
eta3 = ifft(eta_amp,'symmetric');

eta = eta1 + eta2 + eta3;
result(x_index,:) = eta;
x_index = x_index+1;
toc
end
x = 0:0.1:10;
y = (1:1:2000)/20/pi;
y = [fliplr(-y),y];
[xx,yy] = meshgrid(x,y);
result_all = [fliplr(result), result];
surf(yy(2000-100:2000+100,:),xx(2000-100:2000+100,:),result_all(:,2000-100:2000+100)');shading interp;view([0 90]);

