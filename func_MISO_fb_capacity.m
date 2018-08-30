function [ C_fb ] = func_MISO_fb_capacity( B1, B2, h_all, SNR_vector )
% given feedback bits B1 and B2, compute the capacity.
Nt = size(h_all,1);
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
C_fb = mean(C_fb_all,2);
