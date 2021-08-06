clear
addpath('/Users/siwei/Downloads/fsu-program-library-master/barycentric_interp_1d');
x = 0:pi/4:2*pi; 
v = sin(x);
% 将查询点定义为 x 范围内更精细的采样点。

xq = pi/15:pi/15:2*pi;
% 在查询点插入函数并绘制结果。
ni = length(xq);
nd = 4;
figure
% vq1 = interp1(x,v,xq);
vq1 = lagcheby1_interp_1d ( nd, x', v', ni, xq' );
plot(x,v,'o',xq,vq1,':.');
xlim([0 2*pi]);
title('(Default) Linear Interpolation');