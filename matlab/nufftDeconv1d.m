function F = nufftDeconv1d(ftau, tau, M, Mr)
% NUFFTDECONV1D Deconvolution for NUFFT in 1d.
%   F = NUFFTDECONV1D(ftau,tau,Mr) computes the fast Fourier 
%   transform of the function ftau.
%
%   ftau is a periodic fuction of Gauss pulses with spread tau.
%   Mr >  is the number of frequencies on . 
%
    mustBeInteger(M);
    mustBeInteger(Mr);
    if (Mr < M)
        error("M cannot be smaller than Mr");
    end

    % Shifted FFT, computed on oversampled grid
    Ftau = fftshift(fft(ftau,Mr));

    % Shifted wavenumbers
    k = (-M/2:M/2-1)'; 
    
    % Calculate shift in index 
    DM = -M/2-(-Mr/2);
    
    % Index for k on oversampled grid
    I = (1+DM:Mr-DM);

    % Deconvolution
    F = sqrt(pi/tau)*exp(k.^2*tau).*Ftau(I)/Mr;
end
