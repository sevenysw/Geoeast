
function vq1 = barylag_k_matrix(k,l,v,p)
% k : interpolation order
% l : data points
% v : data values
% p : interpolation points
% vq1:interpolation values

ni = length(p);
nx = length(l);

nd = k+1;
wd = ( (-1.0) .^ ( mod ( (0:(nd-1))', 2 ) ) ) ...
.* sin ( (1:2:(2*nd-1))' * pi / ( 2 * nd ) );

M = sparse(zeros(ni,nx));

for j=1:ni
    % xs, vs are the points that used for interpolation 
    % find the nearest nd points.
    kt = k;
    pd = ones(nx-kt,1);
    for s=1:nx-kt
        for i=0:kt
            pd(s) = pd(s) * abs(p(j)-l(s+i));
        end
%         pd(s) = prod(abs(p(j) - l(s:s+kt)));
    end
    [~,sj] = min(pd);
    
    t = wd ./ ( p(j)-l(sj:sj+kt) );
    M(j,sj:sj+kt) = t'/sum(t);
    
end

vq1 = M*v;

end