function [ f ] = fcn_single_nls( u, params )

% Unfold the parameters :
 gam = params.nls.gam;
  mu = params.nls.mu;
xpts = params.geom.xpts;
  dx = params.geom.h;

  % Number of grid - points.
  N = length(xpts); 
  
  d2udx2 = [ u(2) - u(1); diff(u,2); u(N-1) - u(N) ];
  f(1:N) = - 0.5 * d2udx2 / dx^2 + sign(gam) * u.^3 - mu * u;
  
  f = f';
  
end