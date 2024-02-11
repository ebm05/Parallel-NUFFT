function ftau = fgg3d(f, x, Mr, Msp, tau)
% FGG3D Fast Gaussian Gridding (FGG) in 3d.
%   For NUFFT type 1.
%
%   Input:
%       f   : source values at x
%       x   : non-uniform coordinates in [0,2pi)^3
%       Mr  : nr of grid points on the oversampled grid
%       Msp : Spreading width (default=12)
%       tau : Spreading parameter
%
%   Output:
%       ftau : f function convolved with the Gaussian function
%
%   Note: This code is for educational purposes. Loops are not the most
%   efficient way to write MATLAB code. This code is however easily
%   tranlated into other high performance languages such as C or Fortran.
%
%   E Boström, 2024-01-26
if (nargin < 5)
    tau = 
end

N = length(f);
P = 2*Msp;
ftau = zeros(Mr,Mr,Mr);

% Precompute E3
E3(1,1:Msp) = exp(-((pi*(1:Msp)/Mr).^2)/tau); 
E3=[fliplr(E3(1:(Msp-1))),1,E3];

% Loop over all non-uniform points
for n = 1:N
    xn = x(n,:);

    % Find nearest grid point with smallest index values
    m = floor(Mr*xn/(2*pi));
    xi = 2*pi/Mr*m;
    dx = xn-xi;

    % Precompute remaining exponentials before spreading
    E1 = exp(-( dx(1)^2 + dx(2)^2 + dx(3)^2 )/(4*tau));
    E2x = exp(pi*dx(1)/(Mr*tau)).^(-Msp+1:Msp);
    E2y = exp(pi*dx(2)/(Mr*tau)).^(-Msp+1:Msp);
    E2z = exp(pi*dx(3)/(Mr*tau)).^(-Msp+1:Msp);

    % Spreading step
    for ix = 1:P
        for iy = 1:P
            for iz = 1:P                
                % Periodic indices
                jx = mod(m(1)+ix-Msp,Mr)+1;
                jy = mod(m(2)+iy-Msp,Mr)+1;
                jz = mod(m(3)+iz-Msp,Mr)+1;

                % Update convolution result vector
                ftau(jx,jy,jz) = ftau(jx,jy,jz) + ...
                    f(n)*...
                    E1*...
                    E2x(ix)*E2y(iy)*E2z(iz)*...
                    E3(ix)*E3(iy)*E3(iz);
            end
        end
    end
end
end