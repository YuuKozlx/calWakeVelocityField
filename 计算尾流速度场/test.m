
% 用于使用三次样条插值拟合色散关系
clc;clear;close;
% load '../求解特征函数/求解所有的特征函数/eigenfunction.mat'
load eigenfunction.mat


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




