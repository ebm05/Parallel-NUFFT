function F = nufftDeconv3d(ftau, tau, M, R)
% NUFFTDECONV3D Deconvolution for NUFFT in 3d.
%   F = NUFFTDECONV3D(ftau,tau,M,Mr) computes the fast Fourier 
%   transform of the function ftau.
%
%   ftau is a periodic fuction of Gauss pulses with spread tau.
%   Mr >  is the number of frequencies on . 
%
    Mr = M*R;

    % Shifted FFT, computed on oversampled grid
    Ftau = fftshift(fftn(ifftshift(ftau)));
    
    % Calculate shift in index 
    s = -M/2-(-Mr/2);

    % Deconvolution
    F = zeros(M,M,M);
    for i1 = 0:M-1
        for i2 = 0:M-1
            for i3 = 0:M-1
                k1 = i1-M/2;
                k2 = i2-M/2;
                k3 = i3-M/2;                
                F(i1+1,i2+1,i3+1) = sqrt(pi/tau)*exp(tau*k1^2)...
                    *sqrt(pi/tau)*exp(tau*k2^2)...
                    *sqrt(pi/tau)*exp(tau*k3^2)...
                    *Ftau(i1+s,i2+s,i3+s)/(Mr^3);
            end
        end
    end
end
