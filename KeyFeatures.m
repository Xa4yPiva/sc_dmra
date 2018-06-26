function [kf, aNorm] = KeyFeatures(envelope, aThreshold)
%KEYFEATURES Summary of this function goes here
%   Detailed explanation goes here

len = length(envelope);
Nfft = 2 ^ nextpow2(len);

[aNorm, aCN] = SubCoeffsA(envelope);
kf.gammaMax = CoeffGammaMax(aCN, Nfft);
[kf.sigmaAP, kf.sigmaDP] = CoeffSigma(envelope, aNorm, aThreshold);
kf.sigmaAA = CoeffSigmaAA(aCN);
kf.sigmaAF = CoeffSigmaAF(envelope, aNorm, aThreshold);


end

