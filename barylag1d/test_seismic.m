clear
load seismic.mat;

l = (1:8:512)';
nx = length(l);
v = data{2}(:,l);

p = (1:1:l(end))';

figure
imagesc(v,[0 255]);
k = 3;
vq1 = barylag_k_vec(k,l,v,p);
figure
imagesc(vq1,[0 255]);