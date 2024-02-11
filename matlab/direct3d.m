function F = direct3d(f,x,y,z,M1,M2,M3,shift,debug)
%DIRECT3D Direct computation of Fourier coefficients in 3d
%   direct3d(f,x,M) computes the Fourier coefficients directly
%   in the direct naive way.
%   Note that there is also a matlab way to write this routine in an
%   one-line form without loops. The purpose of this routine is educational
%   and for validation.
%
%   Input:
%       f  : source values, f(j), j=1,2,...,N,
%       x  : x coordinates in [0,2pi] or [-pi, pi],
%       y  : y coordinates in [0,2pi] or [-pi, pi],
%       z  : z coordinates in [0,2pi] or [-pi, pi],
%       M1 : number of frequencies in x direction,
%       M2 : number of frequencies in y direction,
%       M3 : number of frequencies in z direction,
%       shift : shifted output? 0 (nonshifted), 1 (shifted, default),
%       debug : 0 (off), 1 (on).
%
%   Output:
%       F : Fourier coefficients F(k1,k2,k3), for
%           k{1,2,3} in 0,...,M{1,2,3} (non-shifted),
%           k{1,2,3} in -M{1,2,3}/2,...,M{1,2,3}/2-1 (shifted).
%
%   Note: This code is for educational purposes. Loops are not the most
%   efficient way to write MATLAB code. This code is however easily
%   tranlated into other high performance languages such as C or Fortran.
%
%   E Boström, 2024-01-26
if (nargin < 9)
    debug = 0;
end
if (nargin < 8)
    shift = 1;
end
if ~isequal(length(f),length(x),length(y),length(z))
    error("The x, y and z arrays must have the same length!");
end

if (mod(M1,2)~=0 || mod(M2,2)~=0 || mod(M3,2)~=0)
    error("M1, M2, M3 must have even integer values.")
end

if (debug)
    fprintf("min(x)=%f\nmax(x)=%f\mmin(y)=%f\nmax(y)=%f\nmin(z)=%f\nmax(z)=%f\n",...
        min(x),max(x),min(y),max(y),min(z),max(z));
end
if(min(x)<-pi || max(x)>=2*pi || min(y)<-pi || max(y)>=2*pi || min(z)<-pi || max(z)>=2*pi)
    error("x,y and z must be in the range [0,2pi] or [-pi,pi]")
end

N = length(f);
F = zeros(M1,M2,M3);
switch shift
    case 0
        if (debug); disp('Using nonshifted'); end
        for k1 = 1:M1
            for k2 = 1:M2
                for k3 = 1:M3
                    for n = 1:N
                        F(k1,k2,k3) = F(k1,k2,k3) +...
                            f(n)*exp(-1i*((k1-1)*x(n)+(k2-1)*y(n)+(k3-1)*z(n)));
                    end
                end
            end
        end
    case 1
        if (debug); disp('Using shifted'); end
        for k1 = -M1/2:(M1/2-1)
            for k2 = -M2/2:(M2/2-1)
                for k3 = -M3/2:(M3/2-1)
                    for n = 1:N
                        F(k1+M1/2+1,k2+M2/2+1,k3+M3/2+1) = F(k1+M1/2+1,k2+M2/2+1,k3+M3/2+1) +...
                            f(n)*exp(-1i*(k1*x(n)+k2*y(n)+k3*z(n)));
                    end
                end
            end
        end
end