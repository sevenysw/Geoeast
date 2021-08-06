clear
l = (0:pi/4:20*pi)'; 
v = [sin(l)];

% 将查询点定义为 x 范围内更精细的采样点。

p = (0:pi/16:20*pi)'+0.05;
% p = pi*2;
% 在查询点插入函数并绘制结果。


k = 3;
tic
vq1 = barylag_k_matrix(k,l,v,p);
toc

figure
plot(l,v,'o');hold on;
plot(p,vq1,':.');

xlim([0 2*pi]);
title('(Default) 3 order Interpolation');