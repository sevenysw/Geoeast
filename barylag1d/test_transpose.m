clear
l = (0:pi/4:2*pi)'; 
v = sin(l);

% 将查询点定义为 x 范围内更精细的采样点。

p = (0:pi/16:2*pi)'+0.05;
p = p(10);
% p = pi*2;
% 在查询点插入函数并绘制结果。

figure
plot(l,v,'o');hold on;
k = 3;
% vq1 = lagcheby1_interp_1d ( length(l), l, v, length(p), p );

vq1 = lagcheby_1point ( k, l, v, length(p), p );
plot(p,vq1,':.');

xlim([0 2*pi]);
title('(Default) 3 order Interpolation');