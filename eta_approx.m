clear;
eta_approximation = pi*sqrt(3)/2*2.^(-2*(1:1:8));
figure,
semilogy(eta_approximation, '-s','linewidth',1)
hold on;
% eta_vector = [0.3634 0.1175 0.03454 0.009497 0.002499 0.0006642 0.0001660 0.00004151];
eta_vector = [0.3634 0.1188 0.03744 0.01154  0.003504 0.001035  0.0002999 0.00008543];
semilogy(eta_vector,'-d','linewidth',1)
grid on;
set(gca, 'fontsize',14)
xlabel('ADC resolution b (bits)')
legend1 = legend('$\frac{\pi \sqrt{3}}{2} 2^{-2b}$', '$\eta_b$');
set(legend1, 'interpreter','latex','fontsize',14);