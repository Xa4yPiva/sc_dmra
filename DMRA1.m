function [ decision ] = DMRA1( keyFeatures, thresholds)
%DMRA1 Summary of this function goes here
%   Detailed explanation goes here

kf = keyFeatures;
thresh = thresholds;

if (kf.gammaMax < thresh.gammaMax)
    if (kf.sigmaAF < thresh.sigmaAF)
        decision = 'FSK2';
    else
        decision = 'FSK4';
    end
else
    if (kf.sigmaAP < thresh.sigmaAP)
        if (kf.sigmaDP < thresh.sigmaDP)
            if (kf.sigmaAA < thresh.sigmaAA)
                decision = 'ASK2';
            else
                decision = 'ASK4';
            end
        else
            decision = 'PSK2';
        end
    else
        decision = 'PSK4';
    end
end

end

