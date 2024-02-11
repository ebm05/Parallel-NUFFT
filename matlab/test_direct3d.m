% Testing the direct 3d sum over an uniformly distributed grid in 3d using 
% all the different shift alternatives:
%   * x,y,z in [0,2pi]  and k1,k2,k3 in {0,...,M-1},
%   * x,y,z in [0,2pi]  and k1,k2,k3 in {-M/2,M/2-1},
%   * x,y,z in [-pi,pi] and k1,k2,k3 in {0,...,M-1},
%   * x,y,z in [-pi,pi] and k1,k2,k3 in {-M/2,M/2-1}.
%
% Note: The x and y coordiates must be shifted for the direct sum, to match
% the way matlab stores three dimensional arrays.
%
% E Boström 2024-01-25

clear
N1 = 16;
N2 = 16;
N3 = 16;
N = N1*N2*N3;
M1 = N1;
M2 = N2;
M3 = N3;
rand_m11 = @(N) 2*rand(N,1)-1;
rel2norm = @(x,y) norm(x-y)/norm(x);

fN = rand_m11(N) + 1i*rand_m11(N);
f = reshape(fN,N1,N2,N3);

% Non-shifted in x, non-shifted in k
x = 2*pi*(0:N1-1)/N1;
y = 2*pi*(0:N2-1)/N2;
z = 2*pi*(0:N3-1)/N3;
[X,Y,Z] = meshgrid(x,y,z);
xN = reshape(X,N,1);
yN = reshape(Y,N,1);
zN = reshape(Z,N,1);
F_dir0 = direct3d(fN,yN,xN,zN,M1,M2,M3,0);
F_FFT0 = fftn(f);
F_dir0 = reshape(F_dir0,N,1);
F_FFT0 = reshape(F_FFT0,N,1);
err00 = rel2norm(F_FFT0,F_dir0);
fprintf("Relative error (2 norm) of direct sum against fftn, non-shifted in x, non-shifted in k: %e\n", err00);

% Non-shifted in x, shifted in k
F_dir01 = direct3d(fN,yN,xN,zN,M1,M2,M3,1);
F_FFT01 = fftshift(fftn(f));
F_dir01 = reshape(F_dir01,N,1);
F_FFT01 = reshape(F_FFT01,N,1);
err01 = rel2norm(F_FFT01,F_dir01);
fprintf("Relative error (2 norm) of direct sum against fftn, non-shifted in x, shifted in k: %e\n", err01);

% Shifted in x, non-shifted in k
x = pi*(2*(0:N1-1)/N1-1);
y = pi*(2*(0:N2-1)/N2-1);
z = pi*(2*(0:N3-1)/N3-1);
[X,Y,Z] = meshgrid(x,y,z);
xN = reshape(X,N,1);
yN = reshape(Y,N,1);
zN = reshape(Z,N,1);
F_dir10 = direct3d(fN,yN,xN,zN,M1,M2,M3,0);
F_FFT10 = fftn(ifftshift(f));
F_dir10 = reshape(F_dir10,N,1);
F_FFT10 = reshape(F_FFT10,N,1);
err10 = rel2norm(F_FFT10,F_dir10);
fprintf("Relative error (2 norm) of direct sum against fftn, shifted in x, shifted in k: %e\n", err10);

% Shifted in x, shifted in k
F_dir11 = direct3d(fN,yN,xN,zN,M1,M2,M3,1);
F_FFT11 = fftshift(fftn(ifftshift(f)));
F_dir11 = reshape(F_dir11,N,1);
F_FFT11 = reshape(F_FFT11,N,1);
err11 = rel2norm(F_FFT11,F_dir11);
fprintf("Relative error (2 norm) of direct sum against fftn, shifted in x, shifted in k: %e\n", err11);

