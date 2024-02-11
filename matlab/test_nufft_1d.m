close all; clear all;

x = [0, pi, 1.8*pi];
f = @(x) cos(x).^2;
Msp = 12;
R = 1;
M = 100;
Mr = R*M;
tau = gaussSpread(M,R,Msp);

xr = linspace(0,2*pi,M);

figure(1)
y = ftau(x,f(x),xr,4,0.01);
plot(xr,y,xr,f(xr));

ytau1 = naiveGridding1d(x,f(x),Mr,Msp,tau);
F1 = nufftDeconv1d(ytau1,tau,M);
figure(2)
plot(xr,ytau1);
title("f_{\tau} for the naive implementation");

ytau2 = fgg1d(x,f(x),Mr,Msp,tau);
F2 = nufftDeconv1d(ytau2,tau,M);
figure(3)
plot(xr,ytau2);
title("f_{\tau} for the FGG implementation")

fprintf("Relative difference between naive and gridding:\nftau:%e\nF:%e.\n",relerr(ytau1,ytau2),relerr(F1,F2));


% nufft in MATLAB
F3 = nufft(f(x),x,xr);


%norm(F1-F2)
%fprintf("Relative difference between naive and gridding:\nFourier coefficients: %e\n",norm(F1-F2)/norm(F1));


% t = [0:300 500.5:700.5];
% S = 2*sin(2*pi*0.02*t) + sin(2*pi*0.1*t);
% X = S + rand(size(t));
% figure(4)
% plot(t,S)
% % figure(2)
% % plot()
% 
% % Create a signal X sampled at unevenly spaced points t. Plot the signal.
% %x = [0:300 500.5:700.5];
% %y = 2*sin(2*pi*0.02*x) + sin(2*pi*0.1*x);
% %X = y + rand(size(t));
% % figure(2)
% % plot(x,y)
% 
% % Compute nufft for the signal
% %N = length(x);
% %M = N;
% %R = 2;
% %Msp = 12;
% %Y = nufft1d_naive(x,y,M,R,Msp)
% 
% 
function y = gtau(x,Msp, tau)
    y = 0*x;
    for l = -Msp:Msp
        y = y + exp(-(x-2*l*pi).^2/(4*tau));
    end   
end

function y = ftau(x,f,xr,Msp,tau)
    N = length(x)
    y = 0*xr;
    for j=1:N
        y = y + f(j)*gtau(x(j)-xr,Msp,tau);
    end
end

function [m, xr] = findNearestBelow(x, Mr)
    m = floor(Mr*x/(2*pi));
    xr = 2*pi*m/Mr;
end
