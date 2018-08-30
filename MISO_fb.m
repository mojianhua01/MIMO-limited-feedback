clear;
close all;
clc;

% use the function func_MISO_fb_capacity to compute the capacity

Nt = 4;
SNR_vector = -20:2:30;
N = 100;
h_all = 1/sqrt(2)* ( randn(Nt, N) + 1j* randn(Nt, N) );

%% capacity with perfect CSIT
for ii = 1:1:length(SNR_vector)
    SNR = 10^(SNR_vector(ii)/10);
    for jj = 1:1:N
        h = h_all(:,jj);
        temp = SNR * norm( h )^2;
        Q = 1 - normcdf( sqrt(temp) );
        C_CSIT_all(ii, jj) = 2* (1 + Q * log_mo(Q) + (1-Q)* log_mo(1-Q) );
    end;
end;
C_CSIT = mean(C_CSIT_all, 2);


%% capacity with limited feedback

B1 = 4;
B2 = 2;

C_fb_0 = func_MISO_fb_capacity(B1, B2, h_all, SNR_vector);

B1 = 4;
B2 = 1;

C_fb_1 = func_MISO_fb_capacity(B1, B2, h_all, SNR_vector);

B1 = 3;
B2 = 1;

C_fb_2 = func_MISO_fb_capacity(B1, B2, h_all, SNR_vector);

B1 = 2;
B2 = 1;

C_fb_3 = func_MISO_fb_capacity(B1, B2, h_all, SNR_vector);


B1 = 1;
B2 = 1;

C_fb_4 = func_MISO_fb_capacity(B1, B2, h_all, SNR_vector);

B1 = 0;
B2 = 0;

C_No_CSIT = func_MISO_fb_capacity(B1, B2, h_all, SNR_vector);

%% plot the figure
figure,
plot(SNR_vector, C_CSIT, 'Linewidth', 1);
hold on;
grid on;
plot(SNR_vector,C_fb_0, 'or--', 'Linewidth',1);
plot(SNR_vector,C_fb_1, 'sr--', 'Linewidth',1);
plot(SNR_vector,C_fb_2, 'dr--', 'Linewidth',1);
plot(SNR_vector,C_fb_3, '^r--', 'Linewidth',1);
plot(SNR_vector,C_fb_4, 'vr--', 'Linewidth',1);
plot(SNR_vector,C_No_CSIT, 'xk-.', 'Linewidth',1.2);
legend('Perfect CSIT','B_1=4, B_2=0', 'B_1=3, B_2=1', 'B_1=2, B_2=2', 'B_1=1, B_2=3','B_1=0, B_2=4', 'No CSIT', 'Location','Best')

set(gca,'FontSize',12);
xlabel('SNR (dB)', 'fontsize', 16)
ylabel('Capacity (bps/Hz)', 'fontsize',16)


figure,
% plot(SNR_vector, C_CSIT, 'Linewidth', 1.2);

plot(SNR_vector,C_CSIT - C_fb_0, 'or-', 'Linewidth',1);
hold on;
grid on;
plot(SNR_vector,C_CSIT - C_fb_1, 'sr-', 'Linewidth',1);
plot(SNR_vector,C_CSIT - C_fb_2, 'dr-', 'Linewidth',1);
plot(SNR_vector,C_CSIT - C_fb_3, '^r-', 'Linewidth',1);
plot(SNR_vector,C_CSIT - C_fb_4, 'vr-', 'Linewidth',1);
plot(SNR_vector,C_CSIT - C_No_CSIT, 'xk-', 'Linewidth',1);
legend('B_1=4, B_2=2', 'B_1=4, B_2=1', 'B_1=3, B_2=1', 'B_1=2, B_2=1','B_1=1, B_2=1', 'No CSIT', 'Location','Best')
set(gca,'FontSize',12);
xlabel('SNR (dB)', 'fontsize', 13)
ylabel('Capacity loss (bps/Hz)', 'fontsize',13)


% %% Bit error rate with perfect CSI
% syms theta;
% f = sqrt(1 - 2/(1 - sin(2*theta)));
% int(f, theta)