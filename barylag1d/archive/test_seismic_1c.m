clear
load seismic.mat;

l = (1:8:512)';
nx = length(l);
v = data{2}(1,l);
% 将查询点定义为 x 范围内更精细的采样点。

p = (1:1:l(end))';
% p = pi*2;
% 在查询点插入函数并绘制结果。

figure
k = 3;
vq1 = barylag_k_vec(k,l,v,p);
plot(l,v,'o',p,vq1,'.-');