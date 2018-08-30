clear;
close all;
M = 10;
B = 10;
h = randn(M,1) + 1j* randn(M,1);
h = h/norm(h);

for i=1:1:2^B
    c(:,i) = randn(M,1) + 1j* randn(M,1);
    c(:,i) = c(:,i)/norm(c(:,i));
end;

inner_prod_square = abs(h'*c).^2;
phase_error = ( angle(h'*c));

figure,
hist(inner_prod_square)
title('Inner product')
figure,
hist(phase_error)
title('Phase')


y_vector = 0:0.001:8;
Q = 1 - normcdf( sqrt(exp(y_vector) ));
Hb = Hb(Q);
figure,
plot(y_vector, Hb, 'linewidth',1.5)
grid on;
hold on;
plot(y_vector, 1 - Hb, 'r--', 'linewidth',1.5)
xlabel('$x$', 'interpreter', 'latex','fontsize',16)
handle = legend('$\mathcal{H}_{\mathrm{b}} \left( Q \left( \sqrt{x} \right) \right)$',...
    '$1-\mathcal{H}_{\mathrm{b}} \left( Q \left( \sqrt{x} \right) \right)$');
set(handle,'Location', 'Best', 'interpreter','latex', 'fontsize', 16)


diff_1_theo = - exp(-y_vector./2)/2./sqrt(2*pi*y_vector).*log2((normcdf(sqrt(y_vector)))./(1 - normcdf(sqrt(y_vector))));
diff_1_simu = diff(Hb,1);
figure,
plot(diff_1_simu/0.01)
title('1st derivative, simulation')
grid on;
figure,
plot(diff_1_theo)
title('1st derivative, theoretical')
grid on;

diff_2_simu =diff(1- Hb,2);
figure,
plot(diff_2_simu/0.001^2)
grid on;
title('2rd derivative, simulation')

clear;
syms x;
Q = 1/2*erfc( x^(1/2) / sqrt(2) );
Hb = -Q.*log2(Q) - (1-Q).*log2(1-Q);
temp = diff(Hb, 2)
figure,
ezplot(temp, x, [0, 10, 0, 10])
% pretty(ans)
% int(Hb, 0, 1)