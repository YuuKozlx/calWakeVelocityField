figure(1)
plot(omega(10000:end),f12(10000:end),"Color","blue")
% 添加0基准线
hold on
plot([0,0],[0,1.2],"Color","black")
hold on;
title('f12')
xlabel('omega')
ylabel('f12')