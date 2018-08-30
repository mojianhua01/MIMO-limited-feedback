clear;
close all;
clc;

% use the function func_MISO_fb_capacity to compute the capacity

Nt = 20;
K = 6
Pt_vector = -20:2:30;
N = 1;
H_all = 1/sqrt(2)* ( randn(K, Nt, N) + 1j* randn(K, Nt, N) );
bit_vector = [1 2 3 4];
eta_vector = [0.3634 0.1175 0.03454 0.009497 0.002499 0.0006642 0.0001660 0.00004151];
bit = 1;

%% capacity with perfect CSIT
for ii = 1:1:length(Pt_vector)
    Pt = 10^(Pt_vector(ii)/10);
    for jj = 1:1:N
        H = H_all(:,:,jj);
        V = pinv(H);
        V_norm = sqrt( sum( abs(V).^2 ) );
        V_normalized = V./V_norm; 
        eta = eta_vector(bit);
        H_V = (H*V_normalized).^2;
        SINR = (1-eta)*Pt/K* diag(H_V)./(eta * Pt/K * diag(H_V) + 1);
        R_all(ii,jj) = mean(log2(1+SINR));
    end;
end;
C_CSIT = mean(C_CSIT_all, 2);


%% capacity with limited feedback

B1 = 4;
B2 = 2;

C_fb_0 = func_MISO_fb_capacity(B1, B2, H_all, Pt_vector);

B1 = 4;
B2 = 1;

C_fb_1 = func_MISO_fb_capacity(B1, B2, H_all, Pt_vector);

B1 = 3;
B2 = 1;

C_fb_2 = func_MISO_fb_capacity(B1, B2, H_all, Pt_vector);

B1 = 2;
B2 = 1;

C_fb_3 = func_MISO_fb_capacity(B1, B2, H_all, Pt_vector);


B1 = 1;
B2 = 1;

C_fb_4 = func_MISO_fb_capacity(B1, B2, H_all, Pt_vector);

B1 = 0;
B2 = 0;

C_No_CSIT = func_MISO_fb_capacity(B1, B2, H_all, Pt_vector);

%% plot the figure
figure,
plot(Pt_vector, C_CSIT, 'Linewidth', 1);
hold on;
grid on;
plot(Pt_vector,C_fb_0, 'or--', 'Linewidth',1);
plot(Pt_vector,C_fb_1, 'sr--', 'Linewidth',1);
plot(Pt_vector,C_fb_2, 'dr--', 'Linewidth',1);
plot(Pt_vector,C_fb_3, '^r--', 'Linewidth',1);
plot(Pt_vector,C_fb_4, 'vr--', 'Linewidth',1);
plot(Pt_vector,C_No_CSIT, 'xk-.', 'Linewidth',1.2);
legend('Perfect CSIT','B_1=4, B_2=0', 'B_1=3, B_2=1', 'B_1=2, B_2=2', 'B_1=1, B_2=3','B_1=0, B_2=4', 'No CSIT', 'Location','Best')

set(gca,'FontSize',12);
xlabel('SNR (dB)', 'fontsize', 16)
ylabel('Capacity (bps/Hz)', 'fontsize',16)


figure,
% plot(SNR_vector, C_CSIT, 'Linewidth', 1.2);

plot(Pt_vector,C_CSIT - C_fb_0, 'or-', 'Linewidth',1);
hold on;
grid on;
plot(Pt_vector,C_CSIT - C_fb_1, 'sr-', 'Linewidth',1);
plot(Pt_vector,C_CSIT - C_fb_2, 'dr-', 'Linewidth',1);
plot(Pt_vector,C_CSIT - C_fb_3, '^r-', 'Linewidth',1);
plot(Pt_vector,C_CSIT - C_fb_4, 'vr-', 'Linewidth',1);
plot(Pt_vector,C_CSIT - C_No_CSIT, 'xk-', 'Linewidth',1);
legend('B_1=4, B_2=2', 'B_1=4, B_2=1', 'B_1=3, B_2=1', 'B_1=2, B_2=1','B_1=1, B_2=1', 'No CSIT', 'Location','Best')
set(gca,'FontSize',12);
xlabel('SNR (dB)', 'fontsize', 13)
ylabel('Capacity loss (bps/Hz)', 'fontsize',13)


% %% Bit error rate with perfect CSI
% syms theta;
% f = sqrt(1 - 2/(1 - sin(2*theta)));
% int(f, theta)