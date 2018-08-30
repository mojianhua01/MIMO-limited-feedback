clc;
clearvars;
% a function mapping from the correlation coefficient of two Gaussian input
% signals to the correlation coefficient of two output signals
Bit_vector = 1:1:8; % [1 2 3 8];
eta_vector = [0.3634 0.1175 0.03454 0.009497 0.002499 0.0006642 0.0001660 0.00004151];
Rho_vector = [0:0.05:0.9, 0.92:0.02:1];
N = 1e5;
for iBit = 1:1:length(Bit_vector)
    bit = Bit_vector(iBit);
    eta = eta_vector(bit);
    for iRho = 1:1:length(Rho_vector)
        rho = Rho_vector(iRho);
%         mu = [0,0];
%         R_yy = [1, rho/2 + 1j*rho/2 ; rho/2 - 1j*rho/2, 1];
        C_yy = [1 rho; rho 1];
        rng default  % For reproducibility
%         y = mvnrnd(mu,R_yy,10000)';
        y = C_yy^(1/2) * randn(2, N);
%         y = C_yy^(1/2)* sqrt(1/2)* ((randn(2,1e4) + 1j* randn(2, 1e4)));
        r = func_quantize(y, bit);
        
        C_yy
        for ii = 1:1:size(r,2)
            C_yy_ii(:,:,ii) = y(:,ii)* y(:,ii)';
        end;
        C_yy_experiment = mean(C_yy_ii,3);
        
        for ii = 1:1:size(r,2)
            C_rr_ii(:,:,ii) = r(:,ii)* r(:,ii)';
        end;
        C_rr = mean(C_rr_ii,3);
        
        Corr = diag(diag(C_rr))^(-1/2) * C_rr * diag(diag(C_rr))^(-1/2);
        
        f(iBit, iRho) = Corr(2,1);
        g(iBit, iRho) = (1-eta)*rho;
        h(iBit, iRho) =  eta * rho^3 + (1-eta)*rho;
        
%         C_rr_theory = (1-eta) * diag(diag((C_yy)))^(1/2) * ...
%     (2/pi * asin( diag(diag((C_yy)))^(-1/2) *...
%     real(C_yy)* diag(diag((C_yy)))^(-1/2) ) + 1j*...
%     2/pi * asin( diag(diag((C_yy)))^(-1/2) *...
%     imag(C_yy)* diag(diag((C_yy)))^(-1/2) ))* ...
%     diag(diag((C_yy)))^(1/2)
        
    end;
end;
% g(:,end) = 1;
figure, p1 = plot(Rho_vector, f(1,:), 'r');
hold on;
plot(Rho_vector, f(2:end,:)', 'r')
p2 = plot(Rho_vector, g(1,:),'--b');
plot(Rho_vector, g(2:end,:)','--b');
% plot(rho_vector, h','-.');
% legend(['1-bit';'2-bit';'3-bit';'4-bit';'5-bit';'6-bit';'7-bit';'8-bit'], 'location','best')
legend([p1; p2], 'f(\phi)', '(1-\eta_b)\phi')
set(gca,'FontSize',12);
xlabel('\phi', 'fontsize', 16)
ylabel('f(\phi)', 'fontsize',16)
grid on;







