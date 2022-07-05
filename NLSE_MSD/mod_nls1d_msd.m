function resid = mod_nls1d_msd( psi, params )
%MOD PARAMS

% Unpack parameters:
npts = params.m;
   dx = params.h;
    V = params.Vstat';
 mu = params.mu;

% Pre-allocate nonlinear residual:
resid = zeros(2*npts,1);

% Decompose the field into real and imaginary parts:
X = psi(1:npts); Y = psi(npts+1:2*npts);

% Compute once the density and common term:
  dens = X.^2 + Y.^2;
comm =  dens + V - mu*ones(npts,1);

% Compute the 1D Laplacians inside the domain (i.e., from j=1 to the Nth
% grid point):
d2Xdx2 = diff(X,2) / dx^2; d2Ydx2 = diff(Y,2) / dx^2;

% Compute the common term (see, the term in the square brackets in Eq.
% (3.4) in Hermano Ricardo's paper)--I am calling it \Omega:
term_l = ( d2Xdx2(1).*X(2)+d2Ydx2(1).*Y(2) ) / dens(2);
term_r = ( d2Xdx2(npts-2).*X(npts-1)+d2Ydx2(npts-2).*Y(npts-1) ) / dens(npts-1);
Omega_l = term_l - 2 * ( dens(2) - dens(1) + V(2) - V(1) );
Omega_r = term_r - 2 * ( dens(npts-1) - dens(npts) + V(npts-1) - V(npts) );

% Compute the 2nd-order derivatives of X and Y and the endpoints:
Xdd_l = Omega_l * X(1); Xdd_r = Omega_r * X(npts);
Ydd_l = Omega_l * Y(1); Ydd_r = Omega_r * Y(npts);

% Update the pre-computed derivatives/Concatenate the vectors:
d2Xdx2 = [ Xdd_l; d2Xdx2; Xdd_r ];
d2Ydx2 = [ Ydd_l; d2Ydx2; Ydd_r ];

% Return the system of nonlinear equations:
resid(1:npts) = -0.5 * d2Xdx2 + comm.*X;
resid(npts+1:2*npts) = -0.5 * d2Ydx2 + comm.*Y;

end