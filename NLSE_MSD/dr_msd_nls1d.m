clearvars; close all; clc; format long;

%% Geometry:
  L = 15; npts = 501;
  x = linspace(-L,L,npts);
dx = x(2)-x(1);

%% Parameters for the 1D GPE:
A = -0.2; B = 1.1;
V = A * sech(B*x/2).^2; V = V';
mu = 0.1;
   params.nls.s   = -1;
   params.nls.a   = 1/2;
% Pack parameters:
params.nls.npts = npts;
     params.nls.x = x;
   params.nls.dx = dx;
    params.nls.V = V;
  params.nls.mu = mu;
%% Newtons
  u0init = [ sqrt(mu) * tanh( sqrt(mu) * x' ) + 0.1*ones(npts,1) ; zeros(npts,1) ];
  
  clc;     
 % Newton:
       u0 = u0init;
       tol = 1d-10;
epsilon = 1d-6;
nmax = 100;
fprintf ('k       ||f(x_k)|| \n')

% For IDRS:
s = 12;               % Size of the shadow space.
tol_lin = 1e-3;    % Tolerance in the linear solver.   
nmax_lin = 1e8; % Maximum number of linear iterations. 
for k=1:nmax
    fx = nls1d_msd( u0, params );
    Jx = jac_nls1d_msd( u0, params );

    % Call IDRS  iterative solver:
% [ p, flag, relres, iter, resvec ] = ...
 %                                  idrs( Jx, -fx, s, tol_lin, nmax_lin ); 
 p = gmres(Jx,-fx,npts,1e-4);
%  p = bicgstab(Jx,-fx,1e-4,100); %***
  %  Jx = dfdjac( @(u_)nls1d_msd(u_,params), u0, epsilon );
  fprintf ('%d     %e     \n',k-1,norm(fx) );
    %p = -Jx \ fx; % Direct Linear Solver.
    norm(p)
  u0 = u0  + p;
%   plot(x,u0(1:npts).^2+u0(npts+1:2*npts).^2);
%   pause;
  if norm(p) < tol*(1+norm(u0)) % norm(fx)<tol%
      break;
  end
end

figure;
u0=u0(1:npts)+1i*u0(npts+1:2*npts);
plot(x,abs(u0).^2,'linewidth',3);
xlabel('$x$','interpreter','latex');
ylabel('$|\psi|^{2}$','interpreter','latex');
set(gca,'fontsize',24,'fontname','times');