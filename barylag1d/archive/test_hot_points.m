clear
addpath('/Users/siwei/Downloads/fsu-program-library-master/barycentric_interp_1d');
l = (0:pi/4:2*pi)'; 
v = sin(l);
nx = length(l);
% ����ѯ�㶨��Ϊ x ��Χ�ڸ���ϸ�Ĳ����㡣

p = (0:pi/16:2*pi)'+0.0;
% p = pi*2;
% �ڲ�ѯ����뺯�������ƽ����
ni = length(p);
k = 3;
figure

vq1 = zeros(ni,1);
for j=1:ni
    % xs, vs are the points that used for interpolation 
    % find the nearest nd points.
    plot(l,v,'o'); hold on;
    kt = k;
    pd = ones(nx-kt,1);
    for s=1:nx-kt
        for i=0:kt
            pd(s) = pd(s) * abs(p(j)-l(s+i));
        end
    end
    [~,sj] = min(pd);
    p1 = sj;
    p2 = sj+kt;
    xs = l(p1:p2);
    vs = v(p1:p2);
    vq1(j) = lagcheby1_interp_1d ( k+1, xs, vs, 1, p(j) );
    xlim([0 2*pi]);
    title('(Default) Linear Interpolation');
    plot(l(p1:p2),v(p1:p2),'x',p(j),vq1(j),':.');
    hold off;
    pause(0.5);
end



