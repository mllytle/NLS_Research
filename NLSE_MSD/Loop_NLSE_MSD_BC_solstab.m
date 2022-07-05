%% Modified 1D NLSE with MSD BC script
% (Newton's Solver for Steady State Sol. and Stability Analysis)
% 6/21/22

clear;close all;clc

%% Params
% params.A = 1/2;
params.a = 1/2;
params.mu = 0.1;
params.gam = 1;
% params.s = -1;

params.m = 501; % number of grid points
params.tol = 1e-8;

A = -0.2; B = 1.1; %***
params.V =@(x) A.*sech(B*x./2).^2; % potential function
params.ns = 1000; % number of newtons steps
plot_flag = 0;

max_leng = 300;

%% Loop Eig Storage
eig_arr = zeros(max_leng-1,2*params.m);
max_eig = zeros(1,max_leng-1);
Lser = 2:max_leng; %***

%% Loop Init
for kk = 0:max_leng-2 %298
L = max_leng-kk;

params.h = 0.5;
xdom = [-L:params.h:L];
params.m = length(xdom);

% params.h = xdom(2)-xdom(1);%(2*L)/(params.m+1); 
% xdom = linspace(-L,L,params.m);

params.xdom = xdom;
params.Vstat = params.V(xdom);

%% NL system of eq. (Get Steady State)
phi = NLSE__MSD_BS(params);
% mu = params.mu;
% phi = [ sqrt(mu) * tanh( sqrt(mu) * params.xdom' ) + ones(params.m,1) ;
%             zeros(params.m,1) ];
%% Stability Analysis
% Construct Stability Matrix
M = stabmat(phi, params);

% Find Eigenvalues
eigs = eig(full(M));
eig_arr(kk+1,1:length(eigs)) = eigs';
max_eig(kk+1) = max(real(eigs));

if plot_flag==2 && 0==mod(kk,10)
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
end % FOR END

figure()
plot(flip(linspace(Lser(1),Lser(end),length(max_eig))),max_eig)

if isreal(1i.*eig_arr)
    disp('No real eigenvalues found.')
else
    disp('Real eigenvalues found.')
end

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
       pcorr = bicgstab(jac_32(phi,params),-NLSE1d_msd(phi,params),bictol,bicsteps); % all mine
%        pcorr = bicgstab(jac_32(phi,params),-mod_nls1d_msd(phi,params),bictol,bicsteps); %stathis funct
%        pcorr = bicgstab(mod_jac_nls1d_msd(phi,params),-NLSE1d_msd(phi,params),bictol,bicsteps); %stathis jac
%        pcorr = bicgstab(mod_jac_nls1d_msd(phi,params),-mod_nls1d_msd(phi,params),bictol,bicsteps); %stathis both
       phi = phi + pcorr;
%        disp(norm(pcorr))
        %disp(norm(NLSE1d_msd(phi,params)))
    
    if norm(pcorr) < tol*(1+norm(phi))
        break;
    end
    
    end

end

function A = stabmat(v,params)
    h = params.h;
    m = params.m;
    mu = params.mu;
    V = params.V;

    A11 = full(spdiags([1i*1/(2*h^2)*ones(m,1) ...
        1i*((mu-1/h^2)*ones(m,1)-2.*abs(conj(v(1:m))).^2-V(v(1:m))) ...
        1i*1/(2*h^2)*ones(m,1)],-1:1,m,m));
    A11(1,1) = 1i*((mu-1/h^2)-2.*abs(v(1)).^2 + v(1)/(2*h^2*v(2)));
    A11(m,m) = 1i*((mu-1/h^2)-2.*abs(v(m)).^2 + v(m)/(2*h^2*v(m-1)));

    A12 = full(spdiags(1i*v(1:m).^2,0,m,m));

    A22 = full(spdiags([-1i*1/(2*h^2)*ones(m,1) ...
        -1i*((mu-1/h^2)*ones(m,1)-2.*abs(conj(v(m+1:2*m))).^2-V(v(m+1:2*m))) ...
        -1i*1/(2*h^2)*ones(m,1)],-1:1,m,m));
    A22(1,1) = -1i*((mu-1/h^2)-2.*abs(v(m+1)).^2 + v(m+1)/(2*h^2*v(m+2)));
    A22(m,m) = -1i*((mu-1/h^2)-2.*abs(v(2*m)).^2 + v(2*m)/(2*h^2*v(2*m-1)));

    A21 = full(spdiags(-1i*v(m+1:2*m).^2,0,m,m));

    A = [A11 A12; A21 A22];
    A(find(isnan(A)))=0;

%     A = [A11 A12; -conj(A12) -A11];
end


