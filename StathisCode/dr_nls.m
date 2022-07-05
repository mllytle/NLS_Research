clear all; close all; format long; clc;

% Parameters of the model:
% %%%%%%%%%%%%%%%%%%%%%%%%
   gam = -1.0;             % Positive / negative values for defocusing / focusing NLS
                            % See the 'fcn_single_nls.m' file for notation/reference. 
    mu =  sign(gam) * 1;    % For the bright soliton, set negative values.
    ic = 'bright';          % 'bright' / 'dark' ( initial guess for Newton ).
solver = 'built-in';           % 'built-in' / 'nsoli' ( nonlinear solver used ).
 
% Geometry:
% %%%%%%%%%%%%%%%%%%%%%%%%
   L = 20;                                            % Lattice size.
  dx = 0.2;                                          % Spatial step-size.
xpts = linspace( - L + dx, L - dx, 2 * L / dx - 1 )'; % Grid - points, x_{i}.                          
   N = length(xpts);                                  % Number of grid points.
   
% Pack the parameters: 
% %%%%%%%%%%%%%%%%%%%%%%%%
% We use a the 'params' structure
% which does the job.
   params.nls.gam = gam;
   params.nls.mu  = mu;
 params.geom.xpts = xpts;
    params.geom.h = dx;
    
% Initial Guess:
% %%%%%%%%%%%%%%%%%%%%%%%%

x0 = 0.0;      % Position of the soliton.
switch(ic)     % Switch between bright and dark soliton 
               % ( note to change the value of mu ).
               
    case('bright')
        
        A = sqrt( -2 * mu );
        
        u0 = A * sech( A * ( xpts - x0 * ones(N,1) ) );
        
    case('dark')
    
        u0 = sqrt(mu) * tanh( sqrt(mu) * ( xpts - x0 * ones(N,1) ) );
        
end

% Nonlinear solvers:
% %%%%%%%%%%%%%%%%%%%%%%%%
switch(solver)
    
    case('built-in') % Built-in "fsolve" nonlinear solver.

 options = optimoptions('fsolve','Display','iter','TolX',1e-13,'TolFun',1e-13);
[ sol, fval, iflag, output ] = fsolve(@(u)fcn_single_nls(u,params), u0, options );

    case('nsoli')    % Nsoli by Keller (Newton-Krylov Method).

% Parameters for nsoli:      
max_iter = 400;
max_iter_linear = 300;
etamax = .9;
lmeth = 1;
restart_limit = 100;
sol_parms = [max_iter,max_iter_linear,etamax,lmeth,restart_limit];
error_flag = 0;
tolerance = 1e-12*[1,1];

% Call nsoli solver:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sol, it_hist, ierr] = nsoli(u0,@(u)fcn_single_nls(u,params),tolerance,sol_parms);

% Please, do check the ierr variable !
end

% Plot the solution against the exact one.
figure(1);
set(gca,'FontSize',16);
plot(xpts,u0,'r','LineWidth',2);
hold on;
plot(xpts,sol,'ob','LineWidth',2);
xlabel('x'); ylabel('u(x)');
legend('Exact','Numerical');
%% LINEAR STABILITY ANALYSIS
% FIRST FORMULATION :
% ===================
       y0 = [sol;zeros(N,1)];
  [ jac ] = jac_nls( y0, params );
[vec,lam] = eig(full(jac));
       ll = diag(lam) * ( -1i );

figure(2); 
set(gca,'FontSize',16);
plot(real(ll),imag(ll),'ob','MarkerSize',12);
xlabel('\lambda_{r}'); ylabel('\lambda_{i}');
%%
% PANOS' ANSATZ :
% ===================
 y0 = [sol;zeros(N,1)];
  [ jac ] = jac_nls( y0, params );
[vec,lam] = eig(full(jac));
       ll = diag(lam);
       
figure(3); 
set(gca,'FontSize',16);
plot(real(ll),imag(ll),'ob','MarkerSize',12);
xlabel('\lambda_{r}'); ylabel('\lambda_{i}');
%% DYNAMICAL EVOLUTION 
ti = 0;
tf = 500;
Nt = 5001;
 t = linspace( ti, tf, Nt );
 
 xpts = params.geom.xpts;
    N = length(xpts);
   y0 = [sol(1:N);zeros(N,1)]; 
      
% Solve the single-NLS using the Dormand and Prince Method :
tic; % Count total execution time.
optiondop = dopset('RelTol',1e-13,'AbsTol',1e-13,'MaxIter',1e10);
[ tout, yout, s ] = dop853( @(t,y0)fcn_single_nls_time(t, y0, params), t, y0, optiondop );
out_stat = ['\n nsteps = ', num2str(s.nsteps),' failed = ', num2str(s.nfailed), ' nfevals = ', num2str(s.nfevals),'\n'];
fprintf(out_stat);
fprintf('\nTime ellapsed = %s\n',datestr(datenum(0,0,0,0,0,toc),'HH:MM:SS'));

% Construct the complex vector :
u = yout(:,1:N) + 1i * yout(:,N+1:2*N);

% Plot the |u|^2 over space-time :
figure(4);
set(gca,'FontSize',16);
imagesc(xpts,tout,abs(u).^2); h = colorbar;
xlabel('x'); ylabel('t');
axis([-L L 0 tout(end)]);
set(h,'FontSize',16);