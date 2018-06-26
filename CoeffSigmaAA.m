function [ sigmaAA ] = CoeffSigmaAA( aCentNorm )
%COEFFSIGMAAA Summary of this function goes here
%   Detailed explanation goes here

N = length(aCentNorm);
sumAA = sum(aCentNorm .^ 2);
sumAbsA = sum(abs(aCentNorm));
sigmaAA = sqrt(sumAA / N - (sumAbsA / N)^2);

% figure(5); plot(abs(aCentNorm)); grid on;

end

