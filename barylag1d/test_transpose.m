clear
l = (0:pi/4:2*pi)'; 
v = sin(l);

% ����ѯ�㶨��Ϊ x ��Χ�ڸ���ϸ�Ĳ����㡣

p = (0:pi/16:2*pi)'+0.05;
p = p(10);
% p = pi*2;
% �ڲ�ѯ����뺯�������ƽ����

figure
plot(l,v,'o');hold on;
k = 3;
% vq1 = lagcheby1_interp_1d ( length(l), l, v, length(p), p );

vq1 = lagcheby_1point ( k, l, v, length(p), p );
plot(p,vq1,':.');

xlim([0 2*pi]);
title('(Default) 3 order Interpolation');