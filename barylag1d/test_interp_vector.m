clear
l = (0:pi/4:2*pi); 
v = [sin(l);cos(l)];

% ����ѯ�㶨��Ϊ x ��Χ�ڸ���ϸ�Ĳ����㡣

p = 0:pi/16:2*pi;
% p = pi*2;
% �ڲ�ѯ����뺯�������ƽ����
i = 2;
figure
plot(l,v,'o');hold on;
k = 3;
vq1 = barylag_k_vec(k,l,v,p);
plot(p,vq1,':.');

xlim([0 2*pi]);
title('(Default) 3 order Interpolation');