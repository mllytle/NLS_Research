%% One Component NLS Mini Excercise
clear
close all

%% Parameters
A = 1;
mu = 1;
xdom = linspace(-20,20,399);
jmax = 399;

%% Problem 2
bright = @(x) A*sech(A*x);
dark = @(x) sqrt(mu)*tanh(sqrt(mu)*x);

% figure(1)
% plot(xdom, bright(xdom))
% grid on
% title('Bright Solution of One-Component NLS Eq. with A = 1','Interpreter','latex')
% xlabel('x','Interpreter','latex')
% ylabel('y','Interpreter','latex')
% 
% figure(2)
% plot(xdom,dark(xdom))
% grid on
% title('Dark Solution of One-Component NLS Eq. with $\mu = 1$','Interpreter','latex')
% xlabel('x','Interpreter','latex')
% ylabel('y','Interpreter','latex')

%% Problem 3
interval = [-20 20];
bv = [0 0];
m = 399;

% Bright Sol
mu = -1/2;
bright_sol = nls_onecomp_fdm(interval, bv, m, 1000, mu, -1);
xdom2 = linspace(-20,20,length(bright_sol));

figure(3)
plot(xdom2, bright_sol)
grid on
title('FDM Bright Solution of One-Component NLS Eq. with $\mu = -1/2$','Interpreter','latex')
xlabel('x','Interpreter','latex')
ylabel('y','Interpreter','latex')

% Dark Sol
% mu = 1;
% dark_sol = nls_onecomp_fdm(interval, bv, m, 200, mu, 1);
% 
% figure(4)
% plot(xdom2, dark_sol)
% grid on
% title('FDM Dark Solution of One-Component NLS Eq. with $\mu = 1$','Interpreter','latex')
% xlabel('x','Interpreter','latex')
% ylabel('y','Interpreter','latex')

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
%         bright = @(x) sech(x);
%         sol = bright(linspace(-20,20,m))';
        
        sol = zeros(m,1);
%         sol(199) = 1;  
% % 
%         l1 = @(x) 0.2*x+1;
%         l2 = @(x) -0.2*x+1;
%         sol(150:199) = l1(linspace(-5,0,50));
%         sol(200:250) = l2(linspace(0,5,51));
%         
%         sol = (0.5*cos(linspace(-20,20,m).*pi./(20)) + 0.5)';

        sol(150:250) = (0.5*cos(linspace(-5,5,101).*pi./(5)) + 0.5)';
        
    %     sol(175:225) = 0.5;
%         sol = ones(m,1);
    elseif gam == 1
%         dark = @(x) tanh(x);
%         sol = dark(linspace(-20,20,m))';
        sol = zeros(m,1);
        sol(1:198) = -1;       
        sol(200:end) = 1;
    end
    
    % Newton's method step
    if gam == -1
        for i = 1:nstep
            sol = sol - ((1/h^2)*jacobian(sol,interval, m, mu, gam))\v_bright(sol,interval,bv,m, mu, gam); %*** changed to dark for now
%             sol = sol - jacobian(sol,interval, m, mu, gam)\v_bright(sol,interval,bv,m, mu, gam);
            plot_test(i,1) = norm(v_bright(sol,interval,bv,m, mu, gam));
        end
        
    elseif gam == 1
        for i = 1:nstep
            sol = sol - ((1/h^2)*jacobian(sol,interval, m, mu, gam))\v_dark(sol,interval,bv,m, mu, gam);
%             sol = sol - jacobian(sol,interval, m, mu, gam)\v_dark(sol,interval,bv,m, mu, gam);
        end
    end
    
    figure(1)
    plot(plot_test)
    xlabel('Number of Iterations')
    ylabel('Value of v_bright (function being minimized)')

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


function y = v_bright(sol,interval,bv,m, mu, gam)
    % Eq. 1.3 function vector for bright sol

    % Init
    y = zeros(m,1);
    h = (interval(2)-interval(1))/(m+1);

    % Lower bv
    y(1) = (sol(2) - 2*sol(1) + bv(1))*1/(h^2) + -2*(gam*sol(1)^2 - mu)*sol(1); %***
%     y(1) = -sol(1) + bv(1);
    
    % Upper bv
    y(m) = (bv(2) - 2*sol(m) + sol(m-1))*1/(h^2) + -2*(gam*sol(m)^2 - mu)*sol(m);
%     y(m) = -sol(m) + bv(2);
    
    % Intermediate bv
    for i=2:m-1
        y(i) = (sol(i+1) - 2*sol(i) + sol(i-1))*1/(h^2) + -2*(gam*sol(i)^2 - mu)*sol(i);
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
