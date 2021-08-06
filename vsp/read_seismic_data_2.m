addpath(genpath('../SegyMAT'));
addpath(genpath('../commonfunction'));

fn = '/Volumes/Yu/BaiduYun/dong_data/anxueyong-1/MSS3D-RL2161-RP6547.sgy';
% [d,H] = ReadSegy(fn);
[d,H] = ReadSegy(fn);


x = [H.SourceX];
y = [H.SourceY];
x = x - min(x(:));
y = y - min(y(:));

d = gain(d,0.004,'agc',0.5,1);
figure;plot(x,y,'.');
figure;imagesc(d);colormap gray;

% D = reshape(d,[size(d,1),112,12]);
% c = 1;
% seishow3D(d,100,-c,c);


% save sea3d.mat D;