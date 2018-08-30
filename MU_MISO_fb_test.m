clear;
close all;
clc;

% use the function func_MISO_fb_capacity to compute the capacity

Nt = 4;
K = 2;
SNR_vector = -20:3:40;
N = 200;
H_all = 1/sqrt(2)* ( randn(K, Nt, N) + 1j* randn(K, Nt, N) );
% bit_vector = [1 2 3 4];
eta_vector = [0.3634 0.1175 0.03454 0.009497 0.002499 0.0006642 0.0001660 0.00004151];
bit = inf;
if bit == inf
    eta = 0;
else
    eta = eta_vector(bit);
end;

B1 = 3;

%% capacity with perfect CSIT, ZF

for jj = 1:1:N
    H = H_all(:,:,jj);
    V = H'*(H*H')^(-1);
    V_norm = sqrt( sum( abs(V).^2 ) );
    V_normalized = bsxfun(@rdivide, V, V_norm);
    H_V = abs(H*V_normalized)^2;
    H_V_CSIT_all(:,:,jj) = H_V;
    for ii = 1:1:length(SNR_vector)
        Pt = 10^(SNR_vector(ii)/10);
        %         SINR = (1-eta)*Pt/K* diag(H_V)./(eta * Pt/K * diag(H_V) + 1);
        SINR = (1-eta)*Pt/K* diag(H_V)./(eta * Pt/K * diag(H_V) + ...
            Pt/K * (sum(H_V,2) - diag(H_V)) + 1);
        R_CSIT_ZF_all(ii,jj) = mean(log2(1+SINR));
    end;
end;
R_CSIT_ZF_mean = real(mean(R_CSIT_ZF_all, 2));

figure, plot(SNR_vector, R_CSIT_ZF_mean,'linewidth', 1)
set(gca, 'fontsize',14)
xlabel('SNR (dB)')
ylabel('Rate (bps/Hz)')
grid on;



% %% capacity with perfect CSIT, MRT
%
% for jj = 1:1:N
%     H = H_all(:,:,jj);
%     V = H';
%     V_norm = sqrt( sum( abs(V).^2 ) );
%     V_normalized = V./V_norm;
%     H_V = (H*V_normalized).^2;
%     for ii = 1:1:length(SNR_vector)
%         Pt = 10^(SNR_vector(ii)/10);
%         %         SINR = (1-eta)*Pt/K* diag(H_V)./(eta * Pt/K * diag(H_V) + 1);
%         SINR = (1-eta)*Pt/K* diag(H_V)./(eta * Pt/K * diag(H_V) + ...
%             Pt/K * (sum(H_V,2) - diag(H_V)) + 1);
%         R_CSIT_MRT_all(ii,jj) = mean(log2(1+SINR));
%     end;
% end;
% R_CSIT_MRT_mean = mean(R_CSIT_MRT_all, 2);

% hold on;, plot(SNR_vector, real(R_CSIT_MRT_mean))
% grid on;


%% capacity with limited feedback, ZF
for jj = 1:1:size(H_all,3)
    H = H_all(:,:,jj);
    
    %% limited feedback
    for kk = 1:1:size(H,1)
        h = H(kk,:);
        %% codebook
        c_all = randn(Nt, 2^B1) + 1j* randn(Nt, 2^B1);
        %         c_all(:,1)  = ctranspose(h);
        c_all_norm = sqrt( sum( abs(c_all).^2 ) );
        c_all_normalized = bsxfun(@rdivide, c_all, c_all_norm);  %normalized
        inner_product = h * c_all_normalized;
        [~, max_inner_product_index] = max ( abs(inner_product) );
        h_hat = c_all_normalized(:, max_inner_product_index);
        H_hat(kk,:) = ctranspose(h_hat);
        %         abs(H(kk,:)*H_hat(kk,:)')/norm(H(kk,:))/norm(H_hat(kk,:))
    end;
    
    %% precoding
    V = H_hat' * (H_hat * H_hat')^(-1);
    V_norm = sqrt( sum( abs(V).^2 ) );
    V_normalized = bsxfun(@rdivide, V, V_norm);
    
    H_V = abs(H*V_normalized).^2;
    H_V_fb_all(:,:,jj) = H_V;
    for ii = 1:1:length(SNR_vector)
        Pt = 10^(SNR_vector(ii)/10);
        %% SINR and Rate
        SINR = (1-eta)*Pt/K* diag(H_V)./(eta * Pt/K * diag(H_V) + ...
            Pt/K * (sum(H_V,2) - diag(H_V)) + 1);
        R_fb_ZF_all(ii,jj) = mean(log2(1+SINR));
    end;
end;

R_fb_ZF_mean = real(mean(R_fb_ZF_all, 2));

hold on;
plot(SNR_vector, R_fb_ZF_mean,'--o','linewidth', 1)
legend1 = legend('Perfect CSIT','Limited feedback');
set(legend1,'fontsize',14,'Location','best')

R_CSIT_ZF_mean - R_fb_ZF_mean


%% for testing
for jj = 1:1:N
    H = H_all(:,:,jj);
    V = randn(Nt, K) + 1j* randn(Nt, K);
    V_norm = sqrt( sum( abs(V).^2 ) );
    V_normalized = bsxfun(@rdivide, V, V_norm);
    H_V = (abs(H*V_normalized)).^2;
    H_V_Random_all(:,:,jj) = H_V;
end;

H_V_CSIT_mean = mean(H_V_CSIT_all,3)
H_V_fb_mean = mean(H_V_fb_all,3)
H_V_Random_mean = mean(H_V_Random_all,3)

%Nt/(Nt-1)*2^(-B1/(Nt-1))
(Nt-K+1)*(1 - 2^(B1)*beta(2^(B1), Nt/(Nt-1)) )
(Nt-K+1)*(1 - 2^(B1)*beta(2^(B1), Nt/(Nt-1)) ) + ...
    (K-1)*2^(B1)*beta(2^(B1), Nt/(Nt-1))
Nt/(Nt-1)*2^(B1)*beta(2^(B1), Nt/(Nt-1))


% plot(SNR_vector, log2(1 + ((1-eta) * 10.^(SNR_vector/10)/K * (Nt-K+1)*(1 - 2^(B1)*beta(2^(B1), Nt/(Nt-1)) ) )...
%     ./(eta * (Nt-K+1)* 10.^(SNR_vector/10)/K* (1 - 2^(B1)*beta(2^(B1), Nt/(Nt-1)) ) + ...
%     (K-1) * 10.^(SNR_vector/10)/K * Nt/(Nt-1)* 2^(B1)* beta(2^B1, Nt/(Nt-1) ) + 1)))

% figure,
% plot(SNR_vector, log2(1 + (1-eta)* 10.^(SNR_vector/10)/K *(Nt-K+1) ./(eta * 10.^(SNR_vector/10)/K * (Nt-K+1)+1)))
% grid on;
% plot(SNR_vector, log2(1/eta)*ones(1, length(SNR_vector)),'--');

% %% capacity with limited feedback, MRT
% for jj = 1:1:size(H_all,3)
%     H = H_all(:,:,jj);
%
%     %% limited feedback
%     for kk = 1:1:size(H,1)
%         h = H(kk,:);
%         inner_product = h * c_all_normalized;
%         [~, max_inner_product_index] = max ( abs(inner_product) );
%         h_hat = c_all_normalized(:, max_inner_product_index);
%         H_hat(kk,:) = ctranspose(h_hat);
%     end;
%
%     %% precoding
%     V = H_hat';
%     V_norm = sqrt( sum( abs(V).^2 ) );
%     V_normalized = V./V_norm;
%     eta = eta_vector(bit);
%     H_V = abs(H*V_normalized).^2;
%     for ii = 1:1:length(SNR_vector)
%         Pt = 10^(SNR_vector(ii)/10);
%         %% SINR and Rate
%         SINR = (1-eta)*Pt/K* diag(H_V)./(eta * Pt/K * diag(H_V) + ...
%             Pt/K * (sum(H_V,2) - diag(H_V)) + 1);
%         R_fb_MRT_all(ii,jj) = mean(log2(1+SINR));
%     end;
% end;
%
% R_fb_MRT_mean = mean(R_fb_MRT_all, 2);

% hold on;
% plot(SNR_vector, real(R_fb_MRT_mean),'--s')


% figure,
% plot(SNR_vector, real(R_CSIT_ZF_mean) - real(R_fb_ZF_mean),'--')
% grid on;
