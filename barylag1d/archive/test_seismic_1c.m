clear
load seismic.mat;

l = (1:8:512)';
nx = length(l);
v = data{2}(1,l);
% ����ѯ�㶨��Ϊ x ��Χ�ڸ���ϸ�Ĳ����㡣

p = (1:1:l(end))';
% p = pi*2;
% �ڲ�ѯ����뺯�������ƽ����

figure
k = 3;
vq1 = barylag_k_vec(k,l,v,p);
plot(l,v,'o',p,vq1,'.-');