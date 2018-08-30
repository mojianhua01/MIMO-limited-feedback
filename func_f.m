function [ R_rr ] = func_f( R_yy, bit)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mu = zeros(1, size(R_yy,1));
sigma = R_yy;
rng default  % For reproducibility
length = 1e4;
% y = mvnrnd(mu,sigma,length)';
y = R_yy^(1/2)* sqrt(1/2)* ((randn(size(R_yy,1),length) + 1j* randn(size(R_yy,1), length)));
[r, eta] = func_quantize(y, bit);
% for ii = 1:1:size(y,2)
%     R_yy_ii(:,:,ii) = y(:,ii)* y(:,ii)';
% end;
% R_yy = mean(R_yy_ii,3)

for ii = 1:1:size(r,2)
    R_rr_ii(:,:,ii) = r(:,ii)* r(:,ii)';
end;
R_rr = mean(R_rr_ii,3);

% R_rr_theory = (1-eta) * diag(diag((R_yy)))^(1/2) * ...
%     (2/pi * asin( diag(diag((R_yy)))^(-1/2) *...
%     real(R_yy)* diag(diag((R_yy)))^(-1/2) ) + 1j*...
%     2/pi * asin( diag(diag((R_yy)))^(-1/2) *...
%     imag(R_yy)* diag(diag((R_yy)))^(-1/2) ))* ...
%     diag(diag((R_yy)))^(1/2)
end

