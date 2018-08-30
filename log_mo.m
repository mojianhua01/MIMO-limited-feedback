function [ output_args ] = log_mo( input_args )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[i1 i2] = size(input_args);
output_args = zeros(i1, i2);
for i=1:1:i1
    for j=1:1:i2
        if input_args(i,j) > 0
            output_args(i,j) = log2 (input_args(i,j));
        else
            outpug_args(i,j) = 0;
        end
    end;
end;

