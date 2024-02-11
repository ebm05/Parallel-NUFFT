function F = naive_nufft3d(f,x,y,z,M1,M2,M3,Msp,R,debug,shift)
%NAIVE_NUFFT3d Naive non-uniform FFT in 3d.
%   naive_nufft3d(f,x,y,z,M1,M2,M3,Msp,R,shift,debug) computes the Fourier 
%   coefficients of non-uniform data. It does gridding over an oversampled
%   domain M1r x M2r x M3r, defined by the R factor (Mjr = R*Mj, j=1,2,3).
%   The spreading Msp is chosen as 12 points in default. For more details 
%   see the paper by Greengaard and Lee (Accelerating the Nonuniform Fast 
%   Fourier Transform SIAM Review).
%
%   Input:
%       f  : source values, f(j), j=1,2,...,N,
%       x  : x coordinates in [0,2pi) or [-pi, pi),
%       y  : y coordinates in [0,2pi) or [-pi, pi),
%       z  : z coordinates in [0,2pi) or [-pi, pi),
%       M1 : number of frequencies in x direction,
%       M2 : number of frequencies in y direction,
%       M3 : number of frequencies in z direction,
%       shift : 0 (non-shifted in xyz, shifted in k, default),
%               1 (shifted in xyz, non-shifted in k),
%               2 (non-shifted Ftau); for testing,
%       debug : 0 (off), 1 (on).
%
%   Output:
%       F : Fourier coefficients F(k1,k2,k3), for
%           k{1,2,3} in -M{1,2,3}/2,...,M{1,2,3}/2-1 (shifted).
%
%   Note: This code is for educational purposes. Loops are not the most
%   efficient way to write MATLAB code. This code is however easily
%   tranlated into other high performance languages such as C or Fortran.
%
%   E Boström, 2024-01-27

max_args = nargin(mfilename);
if (nargin < max_args)
    shift = 0;
end
if (nargin < max_args-1)
    debug = 0;
end
if (nargin < max_args-2)
    R=2;
end
if (nargin < max_args-3)
    Msp = 12;
end
if ~isequal(length(f),length(x),length(y),length(z))
    error("The x, y and z arrays must have the same length!");
end
if ~isequal(M1,M2,M3)
    error("For the current implementation, M1, M2 and M3 must be equal!");
end

N = length(f);
M = M1;
Mr = M*R;
tau = (1/M^2)*(pi*Msp)/(R*(R-0.5));

if (debug) 
    fprintf("(debug) max_args=%d, shift=%d, R=%d, Msp=%d\n",max_args,shift,R,Msp);
    fprintf("(debug) N=%d, M=%d, M1=%d, M2=%d, M3=%d, Mr=%d, tau=%f\n",N,M,M1,M2,M3,Mr,tau); 
end

%% Convolution
a = 2*pi/Mr;
b = 1/(4*tau);
ftau = zeros(Mr,Mr,Mr);
for n=1:N
    xn = x(n);
    yn = y(n);
    zn = z(n);
    
    % Find nearest grid points below
    m1 = floor(xn/a);
    m2 = floor(yn/a);
    m3 = floor(zn/a);
    xi1 = a*m1;
    xi2 = a*m2;
    xi3 = a*m3;

    % Distances
    dx = xn-xi1;
    dy = yn-xi2;
    dz = zn-xi3;

    % Spreading
    for ix = -Msp+1:Msp
        for iy = -Msp+1:Msp
            for iz = -Msp+1:Msp
                jx = mod(m1+ix,Mr);
                jy = mod(m2+iy,Mr);
                jz = mod(m3+iz,Mr);
                if (debug)
                    fprintf("(debug) ix,iy,iz,jx,jy,jz=%d,%d,%d,%d,%d,%d\n",...
                            ix,iy,iz,jx,jy,jz);
                end
                ftau(jx+1,jy+1,jz+1) = ftau(jx+1,jy+1,jz+1) + ...
                    f(n)*exp(-b*((dx-a*ix)^2 + (dy-a*iy)^2 + (dz-a*iz)^2));
            end
        end
    end
end

%% Deconvolution
F = zeros(M1,M2,M3);
sf = 1/(Mr*Mr*Mr);
g = @(k) sqrt(pi/tau)*exp(tau*k.^2);
switch shift
    case 0
        Ftau = fftshift(fftn(ftau))*sf;
    case 1
        Ftau = fftshift(fftn(ifftshift(ftau)))*sf;
    case 2
        Ftau = fftn(ftau)*sf;
end

switch shift
    case 2 % For testing. Non-shifted output of Ftau similar to FFTw

        Ftau = fftshift(Ftau);

        d = Mr/2-M1/2; %shift in k for first element
        i1 = 1;
        for k1 = -M1/2:M1/2-1
            i2 = 1;
            for k2 = -M2/2:M2/2-1
                i3 = 1;
                for k3 = -M3/2:M3/2-1
                    s1 = i1+d;
                    s2 = i2+d;
                    s3 = i3+d;
                    if (debug)
                        fprintf("(debug) k1,k2,k3,s1,s2,s3=%d,%d,%d,%d,%d,%d\n",...
                            k1,k2,k3,s1,s2,s3); 
                    end                    
                    F(i1,i2,i3) = g(k1)*g(k2)*g(k3)*Ftau(s1,s2,s3);
                    i3 = i3 + 1;
                end
                i2 = i2 + 1;
            end
            i1 = i1+1;
        end         

    case {0,1} % shifted output
        d = Mr/2-M1/2; %shift in k for first element
        i1 = 1;
        for k1 = -M1/2:M1/2-1
            i2 = 1;
            for k2 = -M2/2:M2/2-1
                i3 = 1;
                for k3 = -M3/2:M3/2-1
                    s1 = i1+d;
                    s2 = i2+d;
                    s3 = i3+d;
                    if (debug)
                        fprintf("(debug) k1,k2,k3,s1,s2,s3=%d,%d,%d,%d,%d,%d\n",...
                            k1,k2,k3,s1,s2,s3); 
                    end                    
                    F(i1,i2,i3) = g(k1)*g(k2)*g(k3)*Ftau(s1,s2,s3);
                    i3 = i3 + 1;
                end
                i2 = i2 + 1;
            end
            i1 = i1+1;
        end        
end


end