clear
l = (0:pi/4:20*pi)'; 
v = [sin(l)];

% ����ѯ�㶨��Ϊ x ��Χ�ڸ���ϸ�Ĳ����㡣

p = (0:pi/16:20*pi)'+0.05;
% p = pi*2;
% �ڲ�ѯ����뺯�������ƽ����


k = 3;
tic
vq1 = barylag_k_matrix(k,l,v,p);
toc

figure
plot(l,v,'o');hold on;
plot(p,vq1,':.');

xlim([0 2*pi]);
title('(Default) 3 order Interpolation');