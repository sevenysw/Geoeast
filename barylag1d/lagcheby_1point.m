function yi = lagcheby_1point ( nd, xd, yd, ni, xi )

%
  wd = ( (-1.0) .^ ( mod ( (0:(nd-1))', 2 ) ) ) ...
    .* sin ( (1:2:(2*nd-1))' * pi / ( 2 * nd ) );

  t = wd ./ ( xi-xd );
  
  denom = sum(t);

  % the sampled grid and normial grid should not overlapped
  % otherwise yt will contain NaN.
  
  % dot production format
  yi =  t'/denom *yd;
  yt= t/denom*yi; % transpose. 
  
  % matrix multiplication formation
%   yi = diag(1./denom)*t'*yd;
%   yt = t*diag(1./denom)*yi; % transpose. 
  
end
