
% 用于使用三次样条插值拟合色散关系
clc;clear;close;
% load '../求解特征函数/求解所有的特征函数/eigenfunction.mat'
load eigenfunction.mat

a = 0.20;
b = 0.05;
U = 0.5;

% mode 1
wm = eigenfunction(1).w;
Q_mode1 = 4*pi*b^2*a*(sin(wm*a/U) - wm*a/U.*cos(wm*a/U))./(wm*a/U).^3;


k = eigenfunction(1).k;
w1 = eigenfunction(1).w;
w2 = eigenfunction(2).w;
w3 = eigenfunction(3).w;

[fitobject1,~,~] = fit(k,w1,'spline');
[fitobject2,~,~] = fit(k,w2,'spline');
[fitobject3,~,~] = fit(k,w3,'spline');
k = (0.05:0.05:100.05)';
w1_fit = feval(fitobject1,k);
w2_fit = feval(fitobject2,k);
w3_fit = feval(fitobject3,k);
cg1 = diff(w1_fit)./diff(k);
cp1 = w1_fit(1:end-1)./k(1:end-1);

cg2 = diff(w2_fit)./diff(k);
cp2 = w2_fit(1:end-1)./k(1:end-1);

cg3 = diff(w3_fit)./diff(k);
cp3 = w3_fit(1:end-1)./k(1:end-1);
cp = [k(1:end-1),cp1,cp2,cp3];
cg = [k(1:end-1),cg1,cg2,cg3];
w = [w1_fit,w2_fit,w3_fit];
save('w.mat',"w",'-mat')
save('cp.mat',"cp",'-mat')
save('cg.mat',"cg",'-mat')

cp_polyfit = cp(:,2:4);
cg_polyfit = cg(:,2:4);




U = 0.15;

k = k(1:end-1);
for i = 1:10
    phi = 2*pi*i;
    x(:,:,i) = phi*U.*(1-cp_polyfit.*cg_polyfit./U^2)./(k.*(cp_polyfit-cg_polyfit));
    y(:,:,i) = phi*cg_polyfit.*sqrt(1-cp_polyfit.^2/U^2)./(k.*(cp_polyfit-cg_polyfit));
end
y = (y + conj(y))*0.5;
temp = y(:,:,2);
index = sum(temp(:,1)==0)+1;
x = x(index:end,:,:);
y = y(index:end,:,:);

figure(8)
hold on;

for i = 1:10
plot(x(:,1,i),real(y(:,1,i)),'Color','b','LineStyle','-');
plot(x(:,1,i),-real(y(:,1,i)),'Color','b','LineStyle','-');
plot(x(:,2,i),real(y(:,2,i)),'Color','r','LineStyle','-');
plot(x(:,2,i),-real(y(:,2,i)),'Color','r','LineStyle','-');
plot(x(:,3,i),real(y(:,3,i)),'Color','g','LineStyle','-');
plot(x(:,3,i),-real(y(:,3,i)),'Color','g','LineStyle','-');
xlim([0,max(x(:,1,10))])
ylim([-max(x(:,1,3)),max(x(:,1,3))])
end
hold off
figure(9)
 xx = x(:,:,1);
 yy = y(:,:,1);
 plot(xx(:,1),yy(:,1));


 w = w(68:end-1,:);
 ky = sqrt([k(68:end),k(68:end),k(68:end)].^2-(w/U).^2);
