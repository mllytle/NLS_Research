%% Modified 1D NLSE with MSD BC script 
% (Newton's Solver for Steady State Sol. and Stability Analysis)
% 6/20/22

clear;close all;clc

%% Params
% params.A = 1/2;
params.a = 1/2;
params.mu = 0.1; 
    % ^ *** changed to 0.1 from 1 to create stab sol with A = -0.2
params.gam = 1;
% params.s = -1;

params.m = 501; % number of grid points
L = 20;
xdom = linspace(-L,L,params.m);
params.xdom = xdom;
params.h = xdom(2)-xdom(1);%(2*L)/(params.m+1); %*** mod to dr
params.tol = 1e-8;

A = -0.2; B = 1.1;
params.V =@(x) A.*sech(B*x./2).^2; % potential function
params.Vstat = params.V(xdom);
params.ns = 1000; % number of newtons steps
plot_flag = 1;

%% NL system of eq. (Get Steady State)
phi = NLSE__MSD_BS(params);

%% Stability Analysis
% Construct Stability Matrix
M = stabmat(phi, params);

% Find Eigenvalues
eigs = eig(full(M));

%% Plotting
if plot_flag
    figure()
    hold on
    plot (xdom, phi(1:params.m))
    plot (xdom, phi(params.m+1:2*params.m))
    legend ('real','imag')
    title('Solution')
    xlabel('x')
    ylabel('\phi(x)')

    msphi = abs(phi(1:params.m) + 1i.*phi(params.m+1:2*params.m)).^2;
    figure()
    plot (xdom, msphi)
    xlabel('x')
    ylabel('|\phi(x)|^2')
    title('Modulus Squared')

    figure()
    plot(real(eigs),imag(eigs), '.')
    grid on
    xlabel('Real')
    ylabel('Imaginary')
    title('Eigenvalue Spectrum','Interpreter','latex')
end

%% Functions 
function [phi] = NLSE__MSD_BS(params)
    % Get params
    m = params.m;
    mu = params.mu;
    ns = params.ns;
    tol = params.tol;

    % Set Solver Params
    bicsteps = 10000;
    bictol = 1e-7;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ;
    offset = 0;

    % initial guess
%     dark = @(x) sqrt(mu)*tanh(1.2*sqrt(mu)*x)+0.0005;
%     phi_r = dark(linspace(-20,20,m))';
%     phi_i = zeros(1,m)';%+0.1*ones(1,m)'; %check orientation
%     phi = [phi_r;phi_i];
    phi = [ sqrt(mu) * tanh( sqrt(mu) * params.xdom' ) + offset*ones(m,1) ;
            zeros(m,1) ];
    
    % Newtons step
    for ii = 1:ns
       pcorr = bicgstab(jac_32(phi,params),-NLSE1d_msd(phi,params),bictol,bicsteps); % all mine
%        pcorr = bicgstab(jac_32(phi,params),-mod_nls1d_msd(phi,params),bictol,bicsteps); %stathis funct
%        pcorr = bicgstab(mod_jac_nls1d_msd(phi,params),-NLSE1d_msd(phi,params),bictol,bicsteps); %stathis jac
%        pcorr = bicgstab(mod_jac_nls1d_msd(phi,params),-mod_nls1d_msd(phi,params),bictol,bicsteps); %stathis both
       phi = phi + pcorr;
%        disp(norm(pcorr))
        disp(norm(NLSE1d_msd(phi,params)))
    
    if norm(pcorr) < tol*(1+norm(phi))
        break;
    end
    
    end

end

function A = stabmat(v,params)
    h = params.h;
    m = params.m;
    mu = params.mu;

    A11 = full(spdiags([1i*1/(2*h^2)*ones(m,1) 1i*((mu-1/h)*ones(m,1)-2.*abs(v(1:m)).^2) 1i*1/(2*h^2)*ones(m,1)],-1:1,m,m));
    A11(1,1) = 1i*((mu-1/h)-2.*abs(v(1)).^2 + v(1)/(2*h*v(2)));
    A11(m,m) = 1i*((mu-1/h)-2.*abs(v(m)).^2 + v(m)/(2*h*v(m-1)));

    A12 = full(spdiags(1i*v(1:m).^2,0,m,m));

    A22 = full(spdiags([-1i*1/(2*h^2)*ones(m,1) 1i*((mu-1/h)*ones(m,1)-2.*abs(v(m+1:2*m)).^2) -1i*1/(2*h^2)*ones(m,1)],-1:1,m,m));
    A22(1,1) = -1i*((mu-1/h)-2.*abs(v(m+1)).^2 + v(m+1)/(2*h*v(m+2)));
    A22(m,m) = -1i*((mu-1/h)-2.*abs(v(2*m)).^2 + v(2*m)/(2*h*v(2*m-1)));

    A21 = full(spdiags(-1i*v(m+1:2*m).^2,0,m,m));

%     A = [A11 A12; A21 A22];
    A = [A11 A12; -conj(A12) -A11];
end


