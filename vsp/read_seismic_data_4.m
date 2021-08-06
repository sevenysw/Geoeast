addpath(genpath('../SegyMAT'));
addpath(genpath('../commonfunction'));

fn = '/Volumes/Yu/BaiduYun/dong_data/sea-ovt5026.sgy';
% [d,H] = ReadSegy(fn);
[d,H] = ReadSegy(fn,'SkipData',1);


x = [H.GroupX];
y = [H.GroupY];
x = x - min(x(:));
y = y - min(y(:));

figure;plot(x,y,'.');
% d = gain(d,0.004,'agc',0.5,1);
% figure;imagesc(d);colormap gray;

% D = reshape(d,[size(d,1),112,12]);
% c = 1;
% seishow3D(d,100,-c,c);


% save sea3d.mat D;