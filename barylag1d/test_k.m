clear
l = (0:pi/4:2*pi)'; 
v = sin(l);

% ����ѯ�㶨��Ϊ x ��Χ�ڸ���ϸ�Ĳ����㡣

p = 0:pi/16:2*pi;
% p = pi*2;
% �ڲ�ѯ����뺯�������ƽ����

figure
plot(l,v,'o');hold on;
for k=0:5
    vq1 = barylag_k(k,l,v,p);
    plot(p,vq1,':.');

end
xlim([0 2*pi]);
title('(Default) Linear Interpolation');
legend('Original','0','1','2','3','4','5');