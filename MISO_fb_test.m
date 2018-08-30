clear;
% close all;
clc;

% contain the details of the computation of feedback capacity, used for
% testing

Nt = 16;

SNR_vector = -20:2:20;

N = 100;
h_all = 1/sqrt(2)* ( randn(Nt, N) + 1j* randn(Nt, N) );



%% capacity with perfect CSIT
for ii = 1:1:length(SNR_vector)
    SNR = 10^(SNR_vector(ii)/10);
    for jj = 1:1:N
        h = h_all(:,jj);
        temp = SNR * norm( h )^2;
        Q = 1 - normcdf( sqrt(temp) );
        C_all(ii, jj) = 2* (1 + Q * log_mo(Q) + (1-Q)* log_mo(1-Q) );
    end;
end;
C = mean(C_all, 2);

%% capacity with no feedback
B1 = 0;
B2 = 0;

v_all = randn(Nt, 2^B1) + 1j* randn(Nt, 2^B1);
v_all_norm = sqrt( sum( abs(v_all).^2 ) );
v_all_normalized = bsxfun(@rdivide, v_all, v_all_norm);  %normalized

Phi = [0:1:2^B2-1]*( pi/2^(B2+1) ) + pi/2^(B2+2);

for ii = 1:1:length(SNR_vector)
    SNR = 10^(SNR_vector(ii)/10);
    
    for jj = 1:1:size(h_all,2)
        h = h_all(:,jj);
        inner_product = h'*v_all_normalized;
        [max_inner_product, max_inner_product_index] = max ( abs(inner_product) );
        phase_inner_product = angle(h'* v_all_normalized(:, max_inner_product_index) );
        phase_error = min( abs( mod(phase_inner_product, pi/2) - Phi) );
        
        temp1 = SNR * max_inner_product^2 * (1 - sin(2 * phase_error));
        Q1 = 1 - normcdf( sqrt(temp1) );
        temp2 = SNR * max_inner_product^2 * (1 + sin(2 * phase_error));
        Q2 = 1 - normcdf( sqrt(temp2) );
        C_fb_all(ii, jj) = 2 - Hb(Q1) - Hb(Q2);
        
%         temp = SNR * norm(h)^2 * (1 - sin(2 * abs(theta)));
%         Q = 1 - normcdf(sqrt(temp));
%         C_bound_all(ii,jj) = 2 * (1 - Hb(Q));
%         
%         temp3 = SNR * norm(h)^2 * (1 - (sin(2 * abs(theta)))^2);
%         Q3 = 1 - normcdf( sqrt(temp3) );
%         C_test_all(ii,jj) = 2 * (1 - Hb(Q3));
    end;
end;
C_fb_2 = mean(C_fb_all,2);


%% capacity with limited feedback
B1 = 0;
B2 = 16;

v_all = randn(Nt, 2^B1) + 1j* randn(Nt, 2^B1);
v_all_norm = sqrt( sum( abs(v_all).^2 ) );
v_all_normalized = bsxfun(@rdivide, v_all, v_all_norm);  %normalized

Phi = [0:1:2^B2-1]*( pi/2^(B2+1) ) + pi/2^(B2+2);

for ii = 1:1:length(SNR_vector)
    SNR = 10^(SNR_vector(ii)/10);
    
    for jj = 1:1:size(h_all,2)
        h = h_all(:,jj);
        inner_product = h'*v_all_normalized;
        [max_inner_product, max_inner_product_index] = max ( abs(inner_product) );
        phase_inner_product = angle(h'* v_all_normalized(:, max_inner_product_index) );
        phase_error = min( abs( mod(phase_inner_product, pi/2) - Phi) );
        
        temp1 = SNR * max_inner_product^2 * (1 - sin(2 * phase_error));
        Q1 = 1 - normcdf( sqrt(temp1) );
        temp2 = SNR * max_inner_product^2 * (1 + sin(2 * phase_error));
        Q2 = 1 - normcdf( sqrt(temp2) );
        C_fb_all(ii, jj) = 2 - Hb(Q1) - Hb(Q2);
        
%         temp = SNR * norm(h)^2 * (1 - sin(2 * abs(theta)));
%         Q = 1 - normcdf(sqrt(temp));
%         C_bound_all(ii,jj) = 2 * (1 - Hb(Q));
%         
%         temp3 = SNR * norm(h)^2 * (1 - (sin(2 * abs(theta)))^2);
%         Q3 = 1 - normcdf( sqrt(temp3) );
%         C_test_all(ii,jj) = 2 * (1 - Hb(Q3));
    end;
end;
C_fb_1 = mean(C_fb_all,2);
% C_bound = mean(C_bound_all,2);
% C_test = mean(C_test_all,2);




%% plot the figure
figure,
plot(SNR_vector, C, 'Linewidth', 1.2);
hold on;
plot(SNR_vector,C_fb_1, 'or--', 'Linewidth',1.2);
grid on;
plot(SNR_vector,C_fb_2, 'xk-', 'Linewidth',1.2);
% plot(SNR_vector, C_bound, 'k')
% plot(SNR_vector, C_test, '--')
legend('Perfect CSIT',['B_1 =',num2str(B1),', B_2=',num2str(B2)], 'No CSIT', 'Location','Best')
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

% figure,
% bar(theta_all)

%% Bit error rate with perfect CSI
% syms theta;
% f = sqrt(1 - 2/(1 - sin(2*theta)));
% int(f, theta)