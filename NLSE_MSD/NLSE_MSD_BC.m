clear;close all;clc
% Began modifying this on 3/24, but realized no longer works as stable as
% it did prev. so starting new revision (mod prefix) to prevent this from
% happening in the future

%% Params
% params.A = 1/2;
params.a = 1/2;
params.mu = 1;
params.gam = 1;
params.m = 501; % number of grid points
L = 20;
xdom = linspace(-L,L,params.m);
params.xdom = xdom;
params.h = (2*L)/(params.m+1); %*** mod to dr
params.s = -1;
params.V =@(x) 0.2.*sech(1.1*x./2).^2; % potential function

params.ns = 10; % number of newtons steps
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
    s = params.s;
    m = params.m;
    a = params.a;
    mu = params.mu;
    h = params.h;
    ns = params.ns;

    % initial guess
%     bright = @(x) sech(1.2*x) + 0.2;
    dark = @(x) sqrt(mu)*tanh(1.2*sqrt(mu)*x)+0.0005;
    phi_r = dark(linspace(-20,20,m))';
    phi_i = zeros(1,m)';%+0.1*ones(1,m)'; %check orientation
    phi = [phi_r;phi_i];
    
    % Newtons step
    for ii = 1:ns
        pcorr = bicgstab(jac_32(phi,params),f(phi,params),1e-4,100);
%         pcorr = gmres(jac_32(phi,params),-f(phi,params),m,1e-4);
        phi = phi - pcorr;
     %   phi = phi - jac_32(phi,params)\f(phi,params);
    end

end

function y = f(phi,params)
    % Get params *** MAKE DARK
    s = params.s;
    m = params.m;
    a = params.a;
    Vfunct = params.V;
    V = Vfunct(phi);
    mu = params.mu;
    h = params.h;
    phi_r = phi(1:m);
    phi_i = phi(m+1:2*m);

    % Call function to be minimized
    y_r = zeros(m,1);
    y_i = zeros(m,1);
    
    % BCs 
    % BC: real
    N_0 = s*(phi_r(1)^2+phi_i(1)^2) - V(1);
    N_1 = s*(phi_r(2)^2+phi_i(2)^2) - V(2);
    A = (((phi_r(1)+phi_r(3)-2*phi_r(2))/h^2) + phi_r(2) + ...
        ((phi_i(1)+phi_i(3)-2*phi_i(2))/h^2) + phi_i(2)) / (phi_r(2)^2+phi_i(2)^2); %*** is approx right?
    d2y_r1 = (A + N_0 + N_1)*phi_r(1); 
    y_r(1) = d2y_r1-s*(phi_r(1)^2+phi_i(1)^2)*phi_r(1) + (V(1)-mu)*phi_r(1);
    
    N_0 = s*(phi_r(m)^2+phi_i(m)^2) - V(m);
    N_1 = s*(phi_r(m-1)^2+phi_i(m-1)^2) - V(m-1);
    A = (((phi_r(m)+phi_r(m-2)-2*phi_r(m-1))/h^2) + phi_r(m-1) + ...
        ((phi_i(m)+phi_i(m-2)-2*phi_i(m-1))/h^2) + phi_i(m-1)) / (phi_r(m-1)^2+phi_i(m-1)^2); %*** is approx right?
    d2y_r2 = (A + N_0 + N_1)*phi_r(m); 
    y_r(m) = d2y_r2-s*(phi_r(m)^2+phi_i(m)^2)*phi_i(m) + (V(m)-mu)*phi_r(m);

    % BC: imag
    N_0 = s*(phi_r(1)^2+phi_i(1)^2) - V(1);
    N_1 = s*(phi_r(2)^2+phi_i(2)^2) - V(2);
    A = (((phi_r(1)+phi_r(3)-2*phi_r(2))/h^2) + phi_r(2) + ...
        ((phi_i(1)+phi_i(3)-2*phi_i(2))/h^2) + phi_i(2)) / (phi_r(2)^2+phi_i(2)^2);
    d2y_i1 = (A + N_0 + N_1)*phi_i(1);
    y_i(1) = d2y_i1-s*(phi_r(1)^2+phi_i(1)^2)*phi_i(1) + (V(1)-mu)*phi_i(1); %*** check V value
    
    N_0 = s*(phi_r(m)^2+phi_i(m)^2) - V(m);
    N_1 = s*(phi_r(m-1)^2+phi_i(m-1)^2) - V(m-1);
    A = (((phi_r(m)+phi_r(m-2)-2*phi_r(m-1))/h^2) + phi_r(m-1) + ...
        ((phi_i(m)+phi_i(m-2)-2*phi_i(m-1))/h^2) + phi_i(m-1)) / (phi_r(m-1)^2+phi_i(m-1)^2); %*** is approx right?
    d2y_i2 = (A + N_0 + N_1)*phi_i(m); 
    y_i(m) = d2y_i2-s*(phi_r(m)^2+phi_i(m)^2)*phi_i(m) + (V(m)-mu)*phi_i(m);

    
    % Main loop
    for jj = 2:m-1 
        y_r(jj) =  (-a/h^2)*(phi_r(jj-1)+phi_r(jj+1)-2*phi_r(jj)) - ...
            s*(phi_r(jj)^2+phi_i(jj)^2)*phi_r(jj) + (V(jj)-mu)*phi_r(jj);

        y_i(jj) =  (-a/h^2)*(phi_i(jj-1)+phi_i(jj+1)-2*phi_i(jj)) - ...
            s*(phi_r(jj)^2+phi_i(jj)^2)*phi_i(jj) + (V(jj)-mu)*phi_i(jj); %*** check V value
    end

    y = [y_r; y_i];
end
