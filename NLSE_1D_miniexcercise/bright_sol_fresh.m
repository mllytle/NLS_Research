clear
close all

%% Params
A = 1;
mu = -1/2;
gam = -1;
m = 399;

%% Call funct
w1 = nls_1cbs_fdm([-20 20], m, 500, mu, gam);
xdom = linspace(-20,20,length(w1));

%% Fig
bright = @(x) sech(x);

figure(1)
plot(xdom,w1)
hold on
plot(linspace(-20,20,m),bright(linspace(-20,20,m)))
legend('FDM Solution', 'True Solution')

%% Functions
function w = nls_1cbs_fdm(inter, m, ns, mu, gam)
    h = (inter(2) - inter(1))/(m+1);
    
%     w = zeros(m,1);
%     l1 = @(x) 0.2*x+1;
%     l2 = @(x) -0.2*x+1;
%     w(150:199) = l1(linspace(-5,0,50));
%     w(200:250) = l2(linspace(0,5,51));
%     w(150:250) = (0.5*cos(linspace(-5,5,101).*pi./(5)) + 0.5)';
    bright = @(x) sech(1.2*x) + 0.2;
%     bright = @(x) sech(x);
    w = bright(linspace(-20,20,m))';
%     w2 = bright(linspace(-20,20,m))';
%     dbright = @(x) -2.*(gam.*sech(x).^2 - mu).*sech(x);
%     db = dbright(linspace(-20,20,m));

    for ii = 1:ns
        if ii == 1
%             disp(norm(f(w, m, h, mu, gam)-db')/m)
        end
        w = w - jac(w, m, h, mu, gam)\f(w, m, h, mu, gam); %*** 0.5
        residual(ii) = norm(f(w, m, h, mu, gam));%norm(w2-w);
    end
    
    figure(2)
%     semilogy([1:ns],residual)
    plot([1:ns],residual)
    disp(residual(1))
end

function y = f(w, m, h, mu, gam)
    y = zeros(m,1);
    
    y(1) = w(2) + (-2 - 2*gam*(h^2)*w(1)^2 + 2*mu*h^2)*w(1);
    y(m) = w(m-1) + (-2 - 2*gam*(h^2)*w(m)^2 + 2*mu*h^2)*w(m); 
    
    for jj = 2:m-1
       y(jj) =  w(jj-1) + w(jj+1) + (-2 - 2*gam*(h^2)*w(jj)^2 + 2*mu*h^2)*w(jj);
    end
end

function A = jac(w, m, h, mu, gam)
    A = zeros(m,m);
    
    for ll = 1:m
%        A(ll,ll) = (-2 + 2*mu*h^2) - gam*6*(h^2)*w(ll)^2;
       A(ll,ll) = (-2/(h^2) + 2*mu) - gam*6*w(ll)^2;
    end
    
    for kk = 1:m-1
%         A(kk,kk+1)=1;
%         A(kk+1,kk)=1;
        A(kk,kk+1) = 1/(h^2);
        A(kk+1,kk) = 1/(h^2);
    end
end