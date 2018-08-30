clc;
clear all;
close all;
Nr = 2;
Nt = 2;

H_example =[  1.1789 - 0.5500i  -0.6608 + 0.6320i   0.5707 + 1.1739i  -0.7437 - 0.0502i
    -0.1427 - 0.1196i  -0.2934 + 1.6776i   0.4799 - 0.5397i  -0.9118 + 0.1286i
    0.6460 + 1.0750i  -0.2650 - 0.2828i  -1.2834 - 0.5321i  -0.4577 - 0.1047i
    0.6573 + 0.2854i   0.1604 - 0.0104i  -1.0802 - 0.2992i   0.8717 - 0.1929i];
H = H_example(1:Nr, 1:Nt);

H = [1 0; 0 5];

% H = randn(Nr, Nt) + 1j* randn(Nr, Nt);

SNR_dB = [-40:2:20];

SNR = 10.^(SNR_dB/10);

for ite=1:1:length(SNR)
    P = SNR(ite);
    prob = qfunc(sqrt(P/trace((H*H')^(-1))));
    h = -prob*log2(prob) - (1-prob)*log2(1-prob);
    Rate_CI(ite) = 2*Nr*(1 - h);
end;


rho = 1 - 2/pi;
for ite = 1:1:length(SNR)
    P = SNR(ite);
    for ii=1:1:Nr
        beta(ii) = P/Nt* norm(H(ii,:))^2;
        d(ii) = (1 - rho)/(1 + rho* beta(ii));
    end;
    D = diag(d);
    Rate_AQNM(1,ite) = real( log2(det(eye(Nt) + P/Nt * H'*D*H)) );
end;
figure,
plot(SNR_dB, Rate_AQNM)
title('Additive Quantization Noise Model')
xlabel('SNR (dB)')
ylabel('Rate (bps/Hz)')


for ite=1:1:length(SNR)
    P = SNR(ite);
    prob = qfunc(sqrt(P/trace((H*H')^(-1))));
    h = -prob*log2(prob) - (1-prob)*log2(1-prob);
    Rate_CI(ite) = 2*Nr*(1 - h);
end;


for ite = 1:1:length(SNR)
    P = SNR(ite);
    for ii=1:1:Nr
        h_norm(ii) = norm( H(ii,:), 2);
        H_normalized(ii,:) = H(ii,:)/norm( H(ii,:), 2);
    end;
    
    for ii=1:1:Nr
        prob(ii) = qfunc(sqrt(P * h_norm(ii)^2 /trace((H_normalized * H_normalized')^(-1))));
        hb(ii) = -prob(ii)*log2(prob(ii)) - (1-prob(ii))*log2 ( 1-prob(ii) );
        Rate_CI_sub(ii) = 2*(1 - hb(ii));
    end;
    Rate_CI_M(ite) = sum(Rate_CI_sub);
end;
figure,
plot(SNR_dB, Rate_CI_M)
hold on;
plot(SNR_dB, Rate_CI, 'r')
legend('Modified Channel Inversion', 'Channel Inversion',0)
xlabel('SNR (dB)')
ylabel('Rate (bps/Hz)')

%% Upper bound of the capacity
max_sig_val = norm(H, 2);
for ite = 1:1:length(SNR)
    P = SNR(ite);
    pmax = qfunc( sqrt( P* max_sig_val^2/Nr ) );
    if pmax == 0
        temp_max = 0;
    else
        temp_max = -pmax*log2(pmax) - (1-pmax)*log2(1-pmax);
    end;
    Capacity_ub(1,ite) = 2*Nr - 2*Nr*temp_max;
end;
figure,
plot(SNR_dB, Capacity_ub)
title('Upper bound')
xlabel('SNR (dB)')
ylabel('Rate (bps/Hz)')