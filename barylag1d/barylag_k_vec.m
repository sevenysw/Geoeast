
function vq1 = barylag_k_vec(k,l,v,p)
% k : interpolation order
% l : data points
% v : data values
% p : interpolation points
% vq1:interpolation values
v = v';
ni = length(p);
nx = length(l);
nv = size(v,2);
vq1 = zeros(ni,nv);
for j=1:ni
    % xs, vs are the points that used for interpolation 
    % find the nearest nd points.
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
    vs = v(p1:p2,:);
    vq1(j,:) = lagcheby1_interp_1d_vec ( k+1, xs, vs, 1, p(j) );
end
vq1 = vq1';
end