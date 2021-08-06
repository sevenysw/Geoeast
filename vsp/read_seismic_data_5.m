addpath(genpath('../SegyMAT'));
addpath(genpath('../commonfunction'));

fn = '/Volumes/Yu/BaiduYun/dong_data/SEQ027_shot_4000.sgy';
% [d,H] = ReadSegy(fn);
[d,H] = ReadSegy(fn,'SkipData',0);


x = [H.GroupX];
y = [H.GroupY];
x = x - min(x(:));
y = y - min(y(:));
x1 = x(1:480:end);
y1 = y(1:480:end);
figure,plot(x,y,'.');hold on;
plot(x1,y1,'ro');
x1 = x1 - min(x1(:));
x1 = x1/max(x1(:));
figure,plot(x1,ones(length(x1),1),'x');
x2 = linspace(0,1,32);
hold on;plot(x2,ones(length(x2),1)+0.1,'ro');
axis([0 1 0 2])

% % d = gain(d,0.004,'agc',0.5,1);
% figure;imagesc(d);colormap gray;

D = reshape(d,[size(d,1),480,16]);
D = permute(D,[1,3,2]);
c = 100;
seishow3D(D,100,-c,c);

x1 = x1(end:-1:1);
save sea3d_2.mat D x1 x2;