function [ f ] = fcn_single_nls_time( t, y0, params )

% Unfold the parameters :
 gam = params.nls.gam;
xpts = params.geom.xpts;
  dx = params.geom.h;

  % Number of grid - points :
  N = length(xpts); 
  
  % Pre-define the rhs of the equations :
  f = zeros(2*N,1);
  
  % Compose a complex vector :
  u = y0(1:N) + 1i * y0(N+1:2*N);
  
  d2udx2 = [ u(2) - u(1); diff(u,2); u(N-1) - u(N) ];
     rhs = (-1i) * ( - 0.5 * d2udx2 / dx^2 + sign(gam) * abs(u).^2 .* u );
  
  f(1:N) = real(rhs); f(N+1:2*N) = imag(rhs); 
  
end