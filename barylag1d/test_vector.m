clear
x = (0:pi/4:2*pi)'; 
v = [sin(x),cos(x)];
% 将查询点定义为 x 范围内更精细的采样点。

xq = (pi/15:pi/15:2*pi)';
% 在查询点插入函数并绘制结果。
ni = length(xq);
nd = 4;

% vq1 = interp1(x,v,xq);
vq1 = lagcheby1_interp_1d_vec ( nd, x, v, ni, xq );
for i=1:size(v,2)
    figure
    plot(x,v(:,i),'o',xq,vq1(:,i),':.');
    xlim([0 2*pi]);
    title('(Default) Linear Interpolation');
end