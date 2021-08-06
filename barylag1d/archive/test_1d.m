clear
addpath('/Users/siwei/Downloads/fsu-program-library-master/barycentric_interp_1d');
x = 0:pi/4:2*pi; 
v = sin(x);
% ����ѯ�㶨��Ϊ x ��Χ�ڸ���ϸ�Ĳ����㡣

xq = pi/15:pi/15:2*pi;
% �ڲ�ѯ����뺯�������ƽ����
ni = length(xq);
nd = 4;
figure
% vq1 = interp1(x,v,xq);
vq1 = lagcheby1_interp_1d ( nd, x', v', ni, xq' );
plot(x,v,'o',xq,vq1,':.');
xlim([0 2*pi]);
title('(Default) Linear Interpolation');