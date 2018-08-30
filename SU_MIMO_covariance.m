clear;
close all;
clc;

Nr = 4;
Length = 1e5;
Ns = 1;
noise_power = 0;
% H_all = 1/sqrt(2)*(randn(Nr, Nt) + 1j* randn(Nr, Nt));
v = randn(Nr,Ns) + 1j * randn(Nr, Ns);
s = (sqrt(1/2) * (randn(Ns, Length) + 1j * randn(Ns, Length)));
s = real(s); % key point
noise = sqrt(noise_power) * sqrt(1/2) * (randn(Nr, Length) + 1j * randn(Nr, Length));
y = v * s + noise;

for ii = 1:1:Length
    R_yy_ii_separate(:,:,ii) = [real(y(:,ii)); imag(y(:,ii))] * [real(y(:,ii)); imag(y(:,ii))]';
end;
R_yy_separate = mean(R_yy_ii_separate,3)
rho_yy = diag(diag(R_yy_separate))^(-1/2)*R_yy_separate * diag(diag(R_yy_separate))^(-1/2)
for ii = 1:1:Length
    R_yy_ii(:,:,ii) = y(:,ii)* y(:,ii)';
end;
R_yy = mean(R_yy_ii,3)

rho_yy = diag(diag(R_yy))^(-1/2)*R_yy * diag(diag(R_yy))^(-1/2)

bit = 1;

%% stepsize
Lloyd_stepsize_vector = [1.5958 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
stepsize_scale_factor = 1;
eta_vector = [1-2/pi 0.1175 0.03454 0.009497 0.002499 0.0006642 0.0001660 0.00004151];



if bit == +inf
    r = y;
    stepsize = 0;
    eta = 0;
else
    eta = eta_vector(bit);
    for ii=1:1:Nr
        y_ii = y(ii,:);
        stepsize(ii) = stepsize_scale_factor * sqrt(mean(real(y_ii).^2)) * Lloyd_stepsize_vector(bit)+...
            1j * stepsize_scale_factor * sqrt(mean(imag(y_ii).^2))* Lloyd_stepsize_vector(bit);
        %                         stepsize = ( 1.01 * max(abs([real(y); imag(y)])))/(2^(bit-1)); % previous code
        
        %     r_bit_real = min ( max( floor( (real(y)-dither_mean)/stepsize) + 2^(bit - 1), 0), 2^bit-1) ; %[0 2^bit-1]
        %     r_bit_imag = min ( max( floor( (imag(y)-dither_mean)/stepsize )+ 2^(bit - 1), 0), 2^bit-1) ; %[0 2^bit-1]
        %
        %     r_bit = r_bit_real + 1j * r_bit_imag;
        
        r(ii,:) = sign(real(y_ii)) .* ( min( ceil( abs(real(y_ii)) /real(stepsize(ii))) , 2^(bit-1) ) - 1/2 ) * real(stepsize(ii))  + ...
            1j* sign(imag(y_ii)) .* ( min( ceil( abs(imag(y_ii)) /imag(stepsize(ii))) , 2^(bit-1) ) - 1/2 ) * imag(stepsize(ii));
    end;
end;

for ii = 1:1:Length
    R_rr_ii_separate(:,:,ii) = [real(r(:,ii)); imag(r(:,ii))] * [real(r(:,ii)); imag(r(:,ii))]';
end;
R_rr_separate = mean(R_rr_ii_separate,3)
for ii = 1:1:Length
    R_rr_ii(:,:,ii) = r(:,ii)* r(:,ii)';
end;
R_rr = mean(R_rr_ii,3)

rho_rr = diag(diag(R_rr))^(-1/2)*R_rr * diag(diag(R_rr))^(-1/2)

noise_q = r - (1-eta)*y;
for ii = 1:1:Length
    R_qq_ii(:,:,ii) = real(noise_q(:,ii))* real(noise_q(:,ii))' + ...
        1j* imag(noise_q(:,ii))* imag(noise_q(:,ii))';
end;
R_qq_separate = mean(R_qq_ii,3)
for ii = 1:1:Length
    R_qq_ii(:,:,ii) = noise_q(:,ii)* noise_q(:,ii)';
end;
R_qq = mean(R_qq_ii,3)

for ii = 1:1:Length
    R_ry_ii(:,:,ii) = r(:,ii)* y(:,ii)';
end;
R_ry = mean(R_ry_ii,3)

display(R_ry(:)./R_yy(:))

R_yy_separate
rho_yy_separate = diag(diag(R_yy_separate))^(-1/2) * R_yy_separate * diag(diag(R_yy_separate))^(-1/2)

R_rr_separate
rho_rr_separate = diag(diag(R_rr_separate))^(-1/2) * R_rr_separate * diag(diag(R_rr_separate))^(-1/2)

R_rr_computed = (1-eta) * diag(diag(R_yy_separate))^(1/2) * 2 /pi * asin(rho_yy_separate) * diag(diag(R_yy_separate))^(1/2)

% vv = v*v'+noise_power * eye(Nr);
% 
% rho_vv = diag(diag(vv))^(-1/2)*vv * diag(diag(vv))^(-1/2)
% 
% diag(diag(vv))^(1/2) * (2/pi) * ( asin(real(rho_vv)) + 1j *  asin(imag(rho_vv)) ) * diag(diag(vv))^(1/2)
% 
% display((1-eta)^2*R_yy+R_qq)
% 
% R_rr

R_yy = v*v' + noise_power*eye(Nr)

(1-eta)/2 * diag(diag(real(R_yy)))^(1/2) * ...
    2 /pi * asin(diag(diag(real(R_yy)))^(-1/2) *...
    real(R_yy)* diag(diag(real(R_yy)))^(-1/2) ) * ...
    diag(diag(real(R_yy)))^(1/2)

(1-eta) * diag(diag((R_yy)))^(1/2) * ...
    (2/pi * asin( diag(diag((R_yy)))^(-1/2) *...
    real(R_yy)* diag(diag((R_yy)))^(-1/2) ) + 1j*...
    2/pi * asin( diag(diag((R_yy)))^(-1/2) *...
    imag(R_yy)* diag(diag((R_yy)))^(-1/2) ))* ...
    diag(diag((R_yy)))^(1/2)


