function [ sigmaAF ] = CoeffSigmaAF( envelope, aNorm, aThreshold )
%COEFFSIGMAAF Summary of this function goes here
%   Detailed explanation goes here

phDiff = angle(envelope(2 : end) .* conj(envelope(1 : end-1)));
aNorm = aNorm(1:end-1);
phDiffNonWeak = phDiff(aNorm > aThreshold);
phNWCent = phDiffNonWeak - mean(phDiffNonWeak);
C = length(phNWCent);
sumPh2 = sum(phNWCent .^ 2);
sumAbsPh = sum(abs(phNWCent));
sigmaAF = sqrt(sumPh2 / C - (sumAbsPh / C)^2);

% figure(7);
% subplot(2,1,1); plot(phDiff); grid on;
% subplot(2,1,2); plot(abs(phDiff)); grid on;
% 
% figure(8);
% subplot(2,1,1); plot(phDiffNonWeak); grid on;
% subplot(2,1,2); plot(abs(phDiffNonWeak)); grid on;
% 
% figure(9);
% subplot(2,1,1); plot(phNWCent); grid on;
% subplot(2,1,2); plot(abs(phNWCent)); grid on;

end

