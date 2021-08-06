function [t,denom,exact] = precal ( nd, xd, yd, ni, xi )

%*****************************************************************************80
%
%% LAGCHEBY1_INTERP_1D evaluates the Lagrange Chebyshev 1 interpolant.
%
%  Discussion:
%
%    The weight vector WD computed below is only valid if the data points
%    XD are, as expected, the Chebyshev Type 1 points for [-1,+1], or a linearly 
%    mapped version for [A,B].  The XD values may be computed by:
%
%      xd = r8vec_cheby1space ( nd, a, b );
%
%    for instance.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 August 2012
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Jean-Paul Berrut, Lloyd Trefethen,
%    Barycentric Lagrange Interpolation,
%    SIAM Review,
%    Volume 46, Number 3, September 2004, pages 501-517.
%
%  Parameters:
%
%    Input, integer ND, the number of data points.
%    ND must be at least 1.
%
%    Input, real XD(ND,1), the data points.
%
%    Input, real YD(ND,1), the data values.
%
%    Input, integer NI, the number of interpolation points.
%
%    Input, real XI(NI,1), the interpolation points.
%
%    Output, real YI(NI,1), the interpolated values.
%
  wd = ( (-1.0) .^ ( mod ( (0:(nd-1))', 2 ) ) ) ...
    .* sin ( (1:2:(2*nd-1))' * pi / ( 2 * nd ) );
  nv = size(yd,2);
  numer = zeros ( ni, nv );
  denom = zeros ( ni, 1 );
  exact = zeros ( ni, 1 );
  for i=1:ni
      t(i,:) = wd ./ ( xi() - xd );
      denom = sum(t);
      for j = 1 : nd
        exact( xi == xd(j) ) = j;
      end
  end

  yi = t' * yd ./ denom;

  j = find ( exact );
  yi(j,:) = yd(exact(j),:);

  return
end
