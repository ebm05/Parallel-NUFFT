% Testing index shifting
M = 4;
P = 2;
Mr = P*M;

I = -Mr/2:(Mr/2-1);
Is = fftshift(I)

disp("fftshift:")
fftshift(Is)
disp("my fftshift")
myshift(Is,length(Is))

function x = myfftshift(x)
    x = circshift(x,floor(size(x)/2));
end

function n = myshift(k,N)
    n = rem(k+N/2,N);
end