clc;
close all;

[B1, B2] = meshgrid(0:1:20, 0:1:20)

Nt = 10;
figure, 
Z = (1 - 2.^(-B1./(Nt-1))).*(1 - 2.^(-B2) );

figure,
imagesc([0.5:20.5],[0.5:20.5],Z);
colormap(jet);

figure,
contourf(B1, B2, Z, 'ShowText', 'on')
%colormap(jet);
%colorbar;
axis square
xlabel('B_1')
ylabel('B_2')
title('$$f(B_1, B_2) := (1-2^{-\frac{B_1}{N_t-1}}) ( 1 - 2^{-B_2} )$$ when $N_t=10$', 'interpreter', 'latex')
set(gca, 'fontsize', 16)

figure, 
bar3(Z)