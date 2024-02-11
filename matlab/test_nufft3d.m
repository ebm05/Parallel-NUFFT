% Testing the nufft3d routines. 
% Comparison to the direct sum.
% (For testing of the direct sum see test_direct3d.m).
%
% E Boström 2024-01-25

clear
debug = false;
N = 100;
M = 4;
M1 = M; M2 = M; M3 = M;
rand_m11 = @(N) 2*rand(N,1)-1;
rand_01 = @(N) rand(N,1);
rel2norm = @(x,y) norm(x-y)/norm(x);

xN = 2*pi*rand_01(N);
yN = 2*pi*rand_01(N);
zN = 2*pi*rand_01(N);
fN = rand_m11(N) + 1i*rand_m11(N);

F_dir = direct3d(fN,xN,yN,zN,M1,M2,M3);
F_dir = reshape(F_dir,M*M*M,1);
F_nufft_Msp6 = naive_nufft3d(fN,xN,yN,zN,M1,M2,M3,6,2,debug);
F_nufft_Msp6 = reshape(F_nufft_Msp6,M*M*M,1);
err_Msp6 = rel2norm(F_dir,F_nufft_Msp6);
F_nufft_Msp12 = naive_nufft3d(fN,xN,yN,zN,M1,M2,M3,12,2,debug);
F_nufft_Msp12 = reshape(F_nufft_Msp12,M*M*M,1);
err_Msp12 = rel2norm(F_dir,F_nufft_Msp12);
F_nufft_nonshifted_FFT = naive_nufft3d(fN,xN,yN,zN,M1,M2,M3,12,2,false,2);
F_nufft_nonshifted_FFT = reshape(F_nufft_nonshifted_FFT,M*M*M,1);
err_nonshifted_FFT = rel2norm(F_dir,F_nufft_nonshifted_FFT);
fprintf("Relative error (2 norm) of nufft against the direct sum:\n");
fprintf("%e (Msp = 6)\n", err_Msp6);
fprintf("%e (Msp = 12)\n", err_Msp12);
fprintf("%e (Msp = 12, nonshifted FFT, like FFTw)\n", err_nonshifted_FFT);