clear
% close all
clc
%% Params
mu = 1;
m = 399;
interval = [-20 20];
bv = [0 0];
gam = 1;
h = (40)/(m+1);

%% FDM Sol
dark_sol = nls_onecomp_fdm(interval, bv, m, 2, mu, gam);
xdom2 = linspace(-20,20,m);

% dark = @(x) sqrt(mu)*tanh(sqrt(mu)*x);
% dark_sol = dark(xdom2)';

figure(2)
plot(xdom2, dark_sol)
grid on
title('FDM Dark Solution of One-Component NLS Eq. with $\mu = 1$','Interpreter','latex')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')

v = dark_sol;

%% Eig
A12 = full(spdiags(gam.*v.^2,0,m,m));
A11 = full(spdiags([-1/(2*h^2)*ones(m,1) gam.*v.*conj(v) -1/(2*h^2)*ones(m,1)],-1:1,m,m));
A = [A11 A12; -conj(A12) -A11];
% eigs = eig(A)./1i;
load('jac1.mat')
load('sol1.mat')

params.nls.gam = gam;
params.nls.mu = mu;
params.geom.xpts = m;
params.geom.h = (interval(2)-interval(1))/(m+1);
v0 = [sol;zeros(m,1)];
j2 = full(jac_nls2ml( v0, params ));
% j3 = full(spdiags([j2(2,1)*ones(m,1) j2(1,1)*ones(m,1) j2(1,2)*ones(m,1)],-1:1,m,m));

eigs = eig(full(jac))./1i;

figure(4)
plot(eigs, '.')
grid on
xlabel('Real')
ylabel('Imaginary')
title('Eigenvalue Spectrum','Interpreter','latex')

%% Plot modulus squared
t = linspace(0,50,m);
% u = exp(-1i*mu.*t).*v;

for ii = 1:50
    u(ii,:) = exp(-1i*mu.*t(ii)).*v;
end

% figure(2)
% % mesh(xdom,t,umat)
% mesh(xdom,t,u(xdom,t))
% grid on
% xlabel('x')
% ylabel('t')
% zlabel('u')
figure(5)
hold on
for ii = 1:50
    clf
    figure(5)
    plot(xdom2,abs(u(ii,:)).^2)
    ylim([-1 1])
    drawnow
end
xlabel('x')
ylabel('u')
title('Stationary Plot of Mod Squared','Interpreter','latex')
hold off
grid on

%% Functions
function sol = nls_onecomp_fdm(interval, bv, m, nstep, mu, gam)
    % Finite difference method to solve eq. 1.3 using Newton's method
    % Inputs:  boundary conditions and number of grid points
    % OUtputs: Approximate solution to eq. 1.3

    % Make boundary conditions more ascessible
    a = interval(1);
    b = interval(2);

    % Initialize for Newton's
    h = (b-a)/(m+1);
    if gam == -1 % init guesses
        sol = zeros(m,1);
    elseif gam == 1
        dark = @(x) tanh(x);
        sol = dark(linspace(-20,20,m))';
%         sol = zeros(m,1);
%         sol(1:198) = -1;       
%         sol(200:end) = 1;
    end
    
    % Newton's method step

    if gam == -1
        for i = 1:nstep
            sol = sol - ((1/h^2)*jacobian(sol,interval, m, mu, gam))\v_bright(sol,interval,bv,m, mu, gam); %*** changed to dark for now
%             sol = sol - jacobian(sol,interval, m, mu, gam)\v_bright(sol,interval,bv,m, mu, gam);
%             plot_test(i,1) = norm(v_bright(sol,interval,bv,m, mu, gam));
        end
        
    elseif gam == 1
        for i = 1:nstep
            sol = sol - ((1/h^2)*jacobian(sol,interval, m, mu, gam))\v_dark(sol,interval,bv,m, mu, gam);
%             sol = sol - jacobian(sol,interval, m, mu, gam)\v_dark(sol,interval,bv,m, mu, gam);
        end
    end
    
end

function J = jacobian(sol,interval, m, mu, gam)
    % Makes Jacobian matrix

    h = (interval(2)-interval(1))/(m+1);
    J = zeros(m,m);
    
    % main diagonal
    for i = 1:m
        J(i,i) = -2*(1/h^2) - 6*gam*sol(i)^2 + 2*mu; 
    end

    % other diagonals of tridiag
    for i = 1:m-1
%         if gam == -1
%             J(i,i+1) = 1;
%             J(i+1,i) = 1;
%         elseif gam == 1
            J(i,i+1) = 1/h^2;
            J(i+1,i) = 1/h^2;
%         end % ***
    end
end


function y = v_dark(sol,interval,bv,m, mu, gam)
    % Eq. 1.3 function vector for dark sol
    % Init
    y = zeros(m,1);
    h = (interval(2)-interval(1))/(m+1);

    % Lower bv
    y(1) = (sol(2) - 2*sol(1) + sol(2))*1/(h^2) + -2*(gam*sol(1)^2 - mu)*sol(1);
    
    % Upper bv
    y(m) = (sol(m-1) - 2*sol(m) + sol(m-1))*1/(h^2) + -2*(gam*sol(m)^2 - mu)*sol(m);

    % Intermediate bv
    for i=2:m-1
        y(i) = (sol(i+1) - 2*sol(i) + sol(i-1))*1/(h^2) + -2*(gam*sol(i)^2 - mu)*sol(i);
    end
end