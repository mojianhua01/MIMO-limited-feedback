clear;
close all;
clc;


SNR_vector = -10:2:30;

N = 1e3;
h_all = 1/sqrt(2)* ( randn(1, N) + 1j* randn(1, N));


%% capacity with perfect CSIT
for ii = 1:1:length(SNR_vector)
    SNR = 10^(SNR_vector(ii)/10);
    for jj = 1:1:length(h_all)
        h = h_all(jj);
        temp = SNR * norm( h )^2;
        Q = 1 - normcdf( sqrt(temp) );
        C_all(ii, jj) = 2* (1 + Q * log_mo(Q) + (1-Q)* log_mo(1-Q) );
    end;
end;
C = mean(C_all, 2);


%% capacity with limited feedback
B = 2;
theta_all = - pi/2^(B+2) + pi/2^(B+1)* rand(1, N);
for ii = 1:1:length(SNR_vector)
    SNR = 10^(SNR_vector(ii)/10);
    
    for jj = 1:1:length(h_all)
        h = h_all(jj);
       
        theta =  theta_all(jj);
%         theta = pi/8;

        temp1 = SNR * norm(h)^2 * (1 - sin(2 * theta));
        Q1 = 1 - normcdf( sqrt(temp1) );
        temp2 = SNR * norm(h)^2 * (1 + sin(2 * theta));
        Q2 = 1 - normcdf( sqrt(temp2) );
        C_fb_all(ii, jj) = 2 + Q1 * log_mo(Q1) + (1-Q1)* log_mo(1-Q1) + Q2 * log_mo(Q2) + (1-Q2)* log_mo(1-Q2);
        
%         temp = SNR * norm(h)^2 * (1 - sin(2 * abs(theta)));
%         Q = 1 - normcdf(sqrt(temp));
%         C_bound_all(ii,jj) = 2 * (1 - Hb(Q));
%         
%         temp3 = SNR * norm(h)^2 * (1 - (sin(2 * abs(theta)))^2);
%         Q3 = 1 - normcdf( sqrt(temp3) );
%         C_test_all(ii,jj) = 2 * (1 - Hb(Q3));
    end;
end;
C_fb_B_2 = mean(C_fb_all,2);
% C_bound = mean(C_bound_all,2);
% C_test = mean(C_test_all,2);


B = 1;
theta_all = - pi/2^(B+2) + pi/2^(B+1)* rand(1, N);
for ii = 1:1:length(SNR_vector)
    SNR = 10^(SNR_vector(ii)/10);
    
    for jj = 1:1:length(h_all)
        h = h_all(jj);
       
        theta =  theta_all(jj);

        temp1 = SNR * norm(h)^2 * (1 - sin(2 * theta));
        Q1 = 1 - normcdf( sqrt(temp1) );
        temp2 = SNR * norm(h)^2 * (1 + sin(2 * theta));
        Q2 = 1 - normcdf( sqrt(temp2) );
        C_fb_all(ii, jj) = 2 + Q1 * log_mo(Q1) + (1-Q1)* log_mo(1-Q1) + Q2 * log_mo(Q2) + (1-Q2)* log_mo(1-Q2);
        
    end;
end;
C_fb_B_1 = mean(C_fb_all,2);



B = 0;
theta_all = - pi/2^(B+2) + pi/2^(B+1)* rand(1, N);
for ii = 1:1:length(SNR_vector)
    SNR = 10^(SNR_vector(ii)/10);
    
    for jj = 1:1:length(h_all)
        h = h_all(jj);
       
        theta =  theta_all(jj);

        temp1 = SNR * norm(h)^2 * (1 - sin(2 * theta));
        Q1 = 1 - normcdf( sqrt(temp1) );
        temp2 = SNR * norm(h)^2 * (1 + sin(2 * theta));
        Q2 = 1 - normcdf( sqrt(temp2) );
        C_fb_all(ii, jj) = 2 + Q1 * log_mo(Q1) + (1-Q1)* log_mo(1-Q1) + Q2 * log_mo(Q2) + (1-Q2)* log_mo(1-Q2);
        
    end;
end;
C_fb_B_0 = mean(C_fb_all,2);


%% plot the figure
figure,
plot(SNR_vector,C);
hold on;
plot(SNR_vector,C_fb_B_2, 'r^--');
grid on;
plot(SNR_vector, C_fb_B_1, 'rs--')
plot(SNR_vector, C_fb_B_0, 'ro--')
legend('Perfect CSIT','B=2', 'B=1','No CSIT', 'Location','Best')
xlabel('SNR (dB)')
ylabel('Rate (bps/Hz)')
set(gca,'FontSize',12);
% figure,
% plot(10.^(SNR_vector/10),C_fb);
% hold on;
% plot(10.^(SNR_vector/10),C, 'r');
% grid on;
% legend('C_{fb}','C',0)
% xlabel('SNR (linear)')
% ylabel('Rate (bps/Hz)')

figure,
bar(theta_all)

%% Bit error rate with perfect CSI
syms theta;
f = sqrt(1 - 2/(1 - sin(2*theta)));
int(f, theta)