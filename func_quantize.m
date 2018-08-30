function [ r, eta ] = func_quantize( y, bit )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% stepsize
Lloyd_stepsize_vector = [1.5958 0.9957 0.586 0.3352 0.1881 0.1041 0.0569 0.0308];
stepsize_scale_factor = 1;
eta_vector = [1-2/pi 0.1175 0.03454 0.009497 0.002499 0.0006642 0.0001660 0.00004151];

stepsize = zeros(size(y,1), 1);
if bit == +inf
    r = y;
    eta = 0;
else
    eta = eta_vector(bit);
    for ii=1:1:size(y,1)
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

end

