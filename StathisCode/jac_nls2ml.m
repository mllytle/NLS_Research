function [ jac ] = jac_nls2ml( y, params )

% Unfold the parameters :
 gam = params.nls.gam;
  mu = params.nls.mu;
xpts = params.geom.xpts;
  dx = params.geom.h;

  % Number of grid - points.
  N = length(xpts); 

  % Auxiliary stuff:
  one = ones(N,1);
 zero = zeros(N,1);
    u = y(1:N) + 1i * y(N+1:2*N); % Convert to a complex vector.
  
  % Tridiagonal matrix + Free BCs imposed:
     A11 = spdiags([one -2 * one one], -1:1, N, N);
    A11(1,1) = -1; A11(N,N) = -1;
     A11 = -0.5 * A11 / dx^2;
    term = 2 * gam * abs(u).^2  - mu * one;   
     A11 = spdiags(spdiags(A11,0)+term,0,A11);
     A12 = spdiags([zero gam * u.^2 zero],-1:1,N,N);
    
     jac = [ A11, A12; ...
        -conj(A12), -A11  ];
  
end