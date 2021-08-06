addpath(genpath('../SegyMAT'));
addpath(genpath('../commonfunction'));

fn = '/Volumes/Yu/BaiduYun/dong_data/anxueyong-1/MSS3D-Source-30360.sgy';
% [d,H] = ReadSegy(fn);
[d,H] = ReadSegy(fn);


x = [H.GroupX];
y = [H.GroupY];
x = x - min(x(:));
y = y - min(y(:));

figure;plot(x,y,'.');

D = reshape(d,[size(d,1),112,12]);
c = 1;
seishow3D(D,100,-c,c);


save sea3d.mat D;