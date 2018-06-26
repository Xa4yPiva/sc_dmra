function [ decisionNew, meanRat ] = CheckDecision( decision, keyFeatures, keyFeaturesEthalon, modulations )
%CHECKDECISION Summary of this function goes here
%   Detailed explanation goes here

kf = keyFeatures;
if (length(keyFeaturesEthalon) == 1)
    kfEth = keyFeaturesEthalon;
else
    idx = find(modulations == decision);
    kfEth = keyFeaturesEthalon(idx);
end
rat.gammaMax = max(kf.gammaMax, kfEth.gammaMax) / min(kf.gammaMax, kfEth.gammaMax);
rat.sigmaAP = max(kf.sigmaAP, kfEth.sigmaAP) / min(kf.sigmaAP, kfEth.sigmaAP);
rat.sigmaDP = max(kf.sigmaDP, kfEth.sigmaDP) / min(kf.sigmaDP, kfEth.sigmaDP);
rat.sigmaAA = max(kf.sigmaAA, kfEth.sigmaAA) / min(kf.sigmaAA, kfEth.sigmaAA);
rat.sigmaAF = max(kf.sigmaAF, kfEth.sigmaAF) / min(kf.sigmaAF, kfEth.sigmaAF);
ratio = [rat.gammaMax, rat.sigmaAP, rat.sigmaDP, rat.sigmaAA, rat.sigmaAF];
meanRat = mean(ratio(ratio < 1e3));
if meanRat > 10
    decisionNew = 'unknown';
else
    decisionNew = decision;
end

end

