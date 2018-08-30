clear;
close all;
clc;

% use the function func_MISO_fb_capacity to compute the capacity
Nt = 16;
Nr = 4;
SNR_vector = -20:2:20;
N = 1;
H_all = 1/sqrt(2)*(randn(Nr, Nt, N) + 1j* randn(Nr, Nt, N));
% bit_vector = [1 2 3 4];
eta_vector = [0.3634 0.1175 0.03454 0.009497 0.002499 0.0006642 0.0001660 0.00004151];
bit = 2;
eta = eta_vector(bit);

%% codebook at the TX
B1 = 4;
c_all = randn(Nt, 2^B1) + 1j* randn(Nt, 2^B1);
%         c_all(:,1)  = ctranspose(h);
c_all_norm = sqrt( sum( abs(c_all).^2 ) );
c_all_normalized = bsxfun(@rdivide, c_all, c_all_norm);  %normalized

% (sqrt(Nt) + sqrt(Nr))^2

%% complex-valued case
%% capacity with perfect CSIT, MRT
for jj = 1:1:N
    H = H_all(:,:,jj);
    [U S V] = svd(H, 'econ');
    v = V(:,1);
    %     [V D] = eig(H'*H);
    %     maxeigenvalue(jj) = max(diag(D));
    %     v = V(:,end);
    for ii = 1:1:length(SNR_vector)
        Pt = 10^(SNR_vector(ii)/10);
        
        %         for nn = 1:1:100
        %             H_v = H*v;
        %             [U, S, V] = svd((eta * Pt * diag(diag(H_v * H_v')) + eye(Nr))^(-1/2)*H, 'econ');
        %             v = V(:,1);
        %         end;
        
        for kk = 1:1:2^B1
            SINR(kk) = c_all_normalized(:,kk)' * H' *...
                (eta * Pt* diag(diag(H*c_all_normalized(:,kk) *...
                c_all_normalized(:,kk)'*H')) + ...
                eye(Nr))^(-1)*H*c_all_normalized(:,kk);
        end;
        [~, max_index] = max ( abs(SINR) );
        c_best = c_all_normalized(:, max_index);
        c_best_normalized = c_best./norm(c_best);
        
        H_v = H*c_best_normalized;
        SINR1 = (1-eta)*Pt * H_v'*(eta * Pt * diag(diag(H_v * H_v')) + eye(Nr))^(-1)*H_v;
        H_v = H*v;
        SINR2 = (1-eta)*Pt * H_v'*(eta * Pt * diag(diag(H_v * H_v')) + eye(Nr))^(-1)*H_v;
        [SINR, ~] = max([SINR1, SINR2]);
        R_CSIT_MRT_approx_all(ii,jj) = log2(1+SINR);
        
       
        H_v = H*c_best_normalized;
        R_rr = func_f(Pt * H_v * H_v'+ eye(Nr), bit);
        SINR1_lb = (1-eta)*Pt * H_v'*(R_rr - (1-eta)^2 * Pt * H_v * H_v' )^(-1)*H_v;
        
         H_v = H*v;
         R_rr = func_f(Pt * H_v * H_v'+ eye(Nr), bit);
        SINR2_lb = (1-eta)*Pt * H_v'*(R_rr - (1-eta)^2 * Pt * H_v * H_v' )^(-1)*H_v;
        
         [SINR, max_index] = max([SINR1_lb, SINR2_lb]);
        R_CSIT_MRT_lb_all(ii,jj) = log2(1+SINR);
        
        
    end;
end;
R_CSIT_MRT_approx_mean = mean(R_CSIT_MRT_approx_all, 2);
R_CSIT_MRT_lb_mean = mean(R_CSIT_MRT_lb_all, 2);

figure,
plot(SNR_vector, real(R_CSIT_MRT_approx_mean), '-or', 'linewidth',1)
grid on;
hold on;
plot(SNR_vector, real(R_CSIT_MRT_lb_mean), '-sb', 'linewidth',1)

%% capacity with limited feedback, MRT
%% codebook
B1 = 6;
c_all = randn(Nt, 2^B1) + 1j* randn(Nt, 2^B1);
%         c_all(:,1)  = ctranspose(h);
c_all_norm = sqrt( sum( abs(c_all).^2 ) );
c_all_normalized = bsxfun(@rdivide, c_all, c_all_norm);  %normalized

for jj = 1:1:size(H_all,3)
    H = H_all(:,:,jj);
    [U, S, V] = svd(H, 'econ');
    v = V(:,1);
    
    %% precoding
    
    for ii = 1:1:length(SNR_vector)
        Pt = 10^(SNR_vector(ii)/10);
        
        %%iterative algorithm, does not work
        %         for nn = 1:1:100
        %             H_v = H*v;
        %             [U, S, V] = svd((eta * Pt * diag(diag(H_v * H_v')) + eye(Nr))^(-1/2)*H, 'econ');
        %             v = V(:,1);
        %         end;
        
        for kk = 1:1:2^B1
            SINR(kk) = c_all_normalized(:,kk)' * H' *...
                (eta * Pt* diag(diag(H*c_all_normalized(:,kk) *...
                c_all_normalized(:,kk)'*H')) + ...
                eye(Nr))^(-1)*H*c_all_normalized(:,kk);
        end;
        [~, max_index] = max ( abs(SINR) );
        c_best = c_all_normalized(:, max_index);
        c_best_normalized = c_best./norm(c_best);
        
        %%maximize the inner product
        %         inner_product = v' * c_all_normalized;
        %         [~, max_inner_product_index] = max ( abs(inner_product) );
        %         c_best = c_all_normalized(:, max_inner_product_index);
        %         c_best_normalized = c_best./norm(c_best);
        
        %%maximize Hx
        %         inner_product = sum(abs(H * c_all_normalized).^2);
        %         [~, max_inner_product_index] = max ( abs(inner_product) );
        %         c_best = c_all_normalized(:, max_inner_product_index);
        %         c_best_normalized = c_best./norm(c_best);
        
        H_v = H*c_best_normalized;
        
        SINR_approx = (1-eta)*Pt * H_v'*(eta * Pt * diag(diag(H_v * H_v')) + eye(Nr))^(-1)*H_v;
        R_fb_MRT_approx_all(ii,jj) = mean(log2(1+SINR_approx));
        
        R_rr = func_f(Pt * H_v * H_v'+ eye(Nr), bit);
        SINR_lb = (1-eta)*Pt * H_v'*(R_rr - (1-eta)^2 * Pt * H_v * H_v' )^(-1)*H_v;
        R_fb_MRT_lb_all(ii,jj) = mean(log2(1+SINR_lb));
    end;
end;
R_fb_MRT_approx_mean = mean(R_fb_MRT_approx_all, 2);
R_fb_MRT_lb_mean = mean(R_fb_MRT_lb_all, 2);
hold on;
plot(SNR_vector, real(R_fb_MRT_approx_mean),'--or','linewidth',1)
plot(SNR_vector, real(R_fb_MRT_lb_mean),'--sb','linewidth',1)
xlabel('SNR (dB)')
ylabel('Rate (bps/Hz)')
set(gca, 'fontsize',14)
legend({'$\widetilde{R}_{\mathrm{MIMO}}(b, \infty)$', '$\overline{R}_{\mathrm{MIMO}}(b, \infty)$', ...
    '$\widetilde{R}_{\mathrm{MIMO}}(b, B)$', '$\overline{R}_{\mathrm{MIMO}}(b, B)$'}, 'Interpreter', 'latex')

% display(log2(1+(1-eta)/eta*Nr))

% mean(maxeigenvalue)

% for i=1:1:Nr
%     A(:,:,i) = Pt* H(i,:)'*H(i,:);
% end;
% cvx_begin
% variable X(Nt,Nt) symmetric
% variable b(Nr,1)
% maximize sum(b)
% subject to
% for i=1:Nr
% trace( (A(:,:,i) - b(i)* (A(:,:,i)+ eye(Nt)))*X ) >= 0;
% end
% X == semidefinite(Nt);
% cvx_end



% R_yy = Pt*H_v*H_v' + eye(Nr)
% 
% (1-eta)/2 * diag(diag(real(R_yy)))^(1/2) * ...
%     2 /pi * asin(diag(diag(real(R_yy)))^(-1/2) *...
%     real(R_yy)* diag(diag(real(R_yy)))^(-1/2) ) * ...
%     diag(diag(real(R_yy)))^(1/2)
% 
% (1-eta) * diag(diag((R_yy)))^(1/2) * ...
%     (2/pi * asin( diag(diag((R_yy)))^(-1/2) *...
%     real(R_yy)* diag(diag((R_yy)))^(-1/2) ) + 1j*...
%     2/pi * asin( diag(diag((R_yy)))^(-1/2) *...
%     imag(R_yy)* diag(diag((R_yy)))^(-1/2) ))* ...
%     diag(diag((R_yy)))^(1/2)
