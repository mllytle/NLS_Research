%% Modified 1d NLSE script (Newton's Solver)
    % modified for trouble shooting
clear;close all;clc

%% Params
% params.A = 1/2;
params.a = 1/2;
params.mu = 1;
params.gam = 1;
% params.s = -1;

params.m = 501; % number of grid points
L = 20;
xdom = linspace(-L,L,params.m);
params.xdom = xdom;
params.h = xdom(2)-xdom(1);%(2*L)/(params.m+1); %*** mod to dr
params.tol = 1e-8;

A = 0.2; B = 1.1;
params.V =@(x) A.*sech(B*x./2).^2; % potential function
params.Vstat = params.V(xdom);
params.ns = 10000; % number of newtons steps
plot_flag = 1;

%% NL system of eq.
phi = NLSE__MSD_BS(params);

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
end

%% Functions 
function [phi] = NLSE__MSD_BS(params)
    % Get params
    m = params.m;
    mu = params.mu;
    ns = params.ns;
    tol = params.tol;

    % Set Solver Params
    bicsteps = 1000000;
    bictol = 1e-6;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ;
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
%        pcorr = bicgstab(jac_32(phi,params),-NLSE1d_msd(phi,params),bictol,bicsteps); % all mine
%        pcorr = bicgstab(jac_32(phi,params),-mod_nls1d_msd(phi,params),bictol,bicsteps); %stathis funct
%        pcorr = bicgstab(mod_jac_nls1d_msd(phi,params),-NLSE1d_msd(phi,params),bictol,bicsteps); %stathis jac
       pcorr = bicgstab(mod_jac_nls1d_msd(phi,params),-mod_nls1d_msd(phi,params),bictol,bicsteps); %stathis both
       phi = phi + pcorr;
%        disp(norm(pcorr))
        disp(norm(NLSE1d_msd(phi,params)))
    
    if norm(pcorr) < tol*(1+norm(phi))
        break;
    end
    
    end

end


