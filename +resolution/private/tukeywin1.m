function m = tukeywin1(n, nf)
%TUKEYWIN2 Generate a 1-D Tukey window.

i = linspace(0, n, n);

m = (sin(4*pi*i/n)).^2;
m((i > n/nf) & (i < n*(1-1/nf))) = 1;

end
