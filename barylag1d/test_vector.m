clear
x = (0:pi/4:2*pi)'; 
v = [sin(x),cos(x)];
% ����ѯ�㶨��Ϊ x ��Χ�ڸ���ϸ�Ĳ����㡣

xq = (pi/15:pi/15:2*pi)';
% �ڲ�ѯ����뺯�������ƽ����
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