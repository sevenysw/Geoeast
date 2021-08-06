
function vq1 = barylag_k_mat(k,l,v,p)
% k : interpolation order
% l : data points
% v : data values
% p : interpolation points
% vq1:interpolation values
[nt,n2,n3] = size(v);
v = permute(v,[1 3 2]);
v = reshape(v,[nt*n3,n2]);
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
%     vq1(j,:) = lagcheby_1point ( k+1, xs, vs, 1, p(j) );
end
vq1 = vq1';
np  = length(p);
vq1 = reshape(vq1,[nt,n3,np]);
vq1 = permute(vq1,[1,3,2]);

end