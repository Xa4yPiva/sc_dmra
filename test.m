addpath(genpath('../common'));
addpath(genpath('../matlab_utils'));
clear;
% clc;
close all;

fs0 = 12e3;
T = 0.1;
[data, hFilter] = RandomBandLimitedSignal(fs0, T, 20, 300, 3600, 4000, 60, 1, 60, 'uniform');
freqs = 0 : fs0/2-1;
h = freqz(hFilter, length(freqs));
plot(freqs, mag2db(abs(h))); set(gcf, 'color', 'w');
grid on;
xlabel('f, Hz');
ylabel('|H(f)|, dB');
return;
fs = 200e3;
factor = fs / fs0;
[p, q] = rat(factor);
data = resample(data, p, q);
lenData = length(data);
N = 2 ^ nextpow2(lenData);
f = (-N/2 : (N-1)/2) * fs / N;
t = (0 : lenData-1) / fs;

bitsNum = 32;
coeffNorm = 2^bitsNum-1;
dataNorm = uint32((data - min(data)) * coeffNorm);
specData = fftshift(abs(fft(data, N) / (1 * N/2)));

% figure(1);
% subplot(2,1,1); plot(t, data); grid on;
% subplot(2,1,2); plot(f, mag2db(specData)); grid on;

lenFrame = 32;
frame = dataNorm(1 : lenFrame);
% dataBin = de2bi(frame, bitsNum, 'left-msb');
dataBin = randi([0, 1], [lenFrame, bitsNum]);
xBin = double(reshape(dataBin, [1, bitsNum * lenFrame]));
xQ = bi2de(reshape(xBin, [2, (bitsNum * lenFrame) / 2])')';

%% Modulation
fc = fs/4;
M2 = 2;
M4 = 4;
sampleCoeff = 8 * 4;
samplesNum2 = M2 * sampleCoeff;
samplesNum4 = M4 * sampleCoeff;
lenSignal = length(xBin) * samplesNum2;
tSignal = (0 : lenSignal-1) / fs;
offset = exp(1i * 2*pi*fc * tSignal);
Ns = 2 ^ nextpow2(lenSignal);
fsignal = (-Ns/2 : (Ns-1)/2) * fs / Ns;

%--psk--
xPSK2 = SignalOversampleNoFilting(pskmod(xBin, M2), samplesNum2);
xPSK4 = SignalOversampleNoFilting(pskmod(xQ, M4), samplesNum4);
yPSK2 = xPSK2 .* offset;
yPSK4 = xPSK4 .* offset;

% figure(2);
% % subplot(2,2,1); plot(real(yPSK2)); grid on;
% % subplot(2,2,2); plot(real(yPSK4)); grid on;
% snr = 20;
% sPSK2 = awgn(xPSK2, snr);
% sPSK4 = awgn(xPSK4, snr);
% scatterplot(sPSK2); grid on;
% scatterplot(sPSK4); grid on;

%--fsk--
fsep2 = fs / (M2*2);
fsep4 = fs / (M4*2);
xFSK2 = fskmod(xBin, M2, fsep2, samplesNum2, fs);
xFSK4 = fskmod(xQ, M4, fsep4, samplesNum4, fs);
yFSK2 = xFSK2 .* offset;
yFSK4 = xFSK4 .* offset;

xFSK2disc = fskmod(xBin, M2, fsep2, samplesNum2, fs, 'discont');
xFSK4disc = fskmod(xQ, M4, fsep4, samplesNum4, fs, 'discont');
yFSK2disc = xFSK2disc .* offset;
yFSK4disc = xFSK4disc .* offset;

specFSK2 = fftshift(abs(fft(yFSK2, Ns)) / (Ns/2));
specFSK4 = fftshift(abs(fft(yFSK4, Ns)) / (Ns/2));
specFSK2disc = fftshift(abs(fft(yFSK2disc, Ns)) / (Ns/2));
specFSK4disc = fftshift(abs(fft(yFSK4disc, Ns)) / (Ns/2));

% figure(3);
% subplot(2,1,1); plot(fsignal, mag2db(specFSK2)); grid on; ylim([-100, 0]);
% % hold on;    plot(fsignal, mag2db(specFSK2disc));
% subplot(2,1,2); plot(fsignal, mag2db(specFSK4)); grid on; ylim([-100, 0]);
% % hold on;    plot(fsignal, mag2db(specFSK4disc));

n = 20;
% figure(4);
% subplot(2,2,1); plot(xBin(1:n)); grid on;
% subplot(2,2,3); plot(xQ(1:n)); grid on;
% subplot(2,2,1); plot(real(yFSK2disc(1 : n * samplesNum2))); grid on;
% subplot(2,2,2); plot(real(yFSK4disc(1 : n * samplesNum4))); grid on;
% subplot(2,2,3); plot(real(yFSK2(1 : n * samplesNum2))); grid on;
% subplot(2,2,4); plot(real(yFSK4(1 : n * samplesNum4))); grid on;

%--ask--
yASK2 = ASKmod(xBin, M2, samplesNum2, fc, fs, [0.2, 1]);
yASK4 = ASKmod(xQ, M4, samplesNum4, fc, fs, [0.25, 0.5, 0.75, 1]);
specASK2 = fftshift(abs(fft(yASK2, Ns)) / (Ns/2));
specASK4 = fftshift(abs(fft(yASK4, Ns)) / (Ns/2));

% figure(5);
% subplot(2,1,1); plot(yASK2(1 : n * nsamp2)); grid on;
% subplot(2,1,2); plot(yASK4(1 : n * nsamp4)); grid on;
% figure(6);
% subplot(2,1,1); plot(fsignal, mag2db(specASK2)); grid on; %ylim([-100, 0]);
% subplot(2,1,2); plot(fsignal, mag2db(specASK4)); grid on; %ylim([-100, 0]);

snr = 7;
signal = yFSK4;
envel = hilbert(real(signal)) .* conj(offset);
envel = awgn(envel(1:4*4096), snr);
% % envel = hilbert(awgn(real(signal), snr)) .* conj(offset);

% figure(8);
% plot(angle(envel(2 : end) .* conj(envel(1:end-1))));
% grid on;

% figure(6);
% % subplot(2,1,1); plot(real(envel(1 : n * samplesNum2))); grid on;
% subplot(2,1,1); plot(real(envel)); grid on;
% subplot(2,1,2); plot(fsignal, mag2db(fftshift(abs(fft(envel, Ns)) / (Ns/2)))); grid on;

% phDiff = angle(envel(2 : end) .* conj(envel(1 : end-1)));
% figure(7);
% plot(phDiff(1 : n * samplesNum2) / pi * fs/2);
% grid on;

modulations = ["ASK2", "ASK4", "PSK2", "PSK4", "FSK2", "FSK4"];

%% DMRA
% Dirty hardcode!
thresholds.ampl     = 1;
thresholds.gammaMax = 0.8;        % 4
thresholds.sigmaAP  = pi/5.5;   % pi/5.5
thresholds.sigmaDP  = pi/3;     % pi/6.5 - pi/2.5
thresholds.sigmaAA  = 0.2;     % 0.25
thresholds.sigmaAF  = 0.4;      % 0.4

thresholds

% kfEth = KeyFeatures(envel, thresholds.ampl);
% envel = awgn(envel, snr);
kf = KeyFeatures(envel, thresholds.ampl)
dmra1 = DMRA1(kf, thresholds);
fprintf("dmra1 = " + dmra1 + "\n");
% [decision, meanRat] = CheckDecision(dmra1, kf, kfEth);

% rat.gammaMax = max(kf.gammaMax, kfEth.gammaMax) / min(kf.gammaMax, kfEth.gammaMax);
% rat.sigmaAP = max(kf.sigmaAP, kfEth.sigmaAP) / min(kf.sigmaAP, kfEth.sigmaAP);
% rat.sigmaDP = max(kf.sigmaDP, kfEth.sigmaDP) / min(kf.sigmaDP, kfEth.sigmaDP);
% rat.sigmaAA = max(kf.sigmaAA, kfEth.sigmaAA) / min(kf.sigmaAA, kfEth.sigmaAA);
% rat.sigmaAF = max(kf.sigmaAF, kfEth.sigmaAF) / min(kf.sigmaAF, kfEth.sigmaAF);
% rat

%% Signals, Decisions
signals = [yASK2; yASK4; yPSK2; yPSK4; yFSK2; yFSK4];
decRight = ["ASK2", "ASK4", "PSK2", "PSK4", "FSK2", "FSK4"];
sigsNum = min(size(signals));
mType = ['x', '*', 's', '^', 'o', 'd-'];

lenF = 2^15;
fNum = floor(max(size(signals)) / lenF);


%% Probability of right decision
% envelopes = zeros(sigsNum, lenSignal);
% for i = 1 : length(modulations)
%     envelopes(i, :) = hilbert(real(signals(i, :))) .* conj(offset);
%     kfEthalon(i) = KeyFeatures(envelopes(i, :), thresholds.ampl);
% end
% SNR = (-2 : 0.5 : 13);
% expNum = 100;
% 
% lenSNR = length(SNR);
% pRight = zeros(sigsNum, lenSNR);
% cyclesNum = sigsNum * lenSNR;
% iteration = 0;
% h = waitbar(0, 'Computing...');
% tic
% for k = 1 : sigsNum
%     for i = 1 : lenSNR
%         decRightNum = 0;
%         for j = 1 : expNum
%             pos = floor(rand() * (fNum-1));
%             env = awgn(envelopes(k, pos*lenF+1 : (pos+1)*lenF), SNR(i), 'measured');
%             kf = KeyFeatures(env, thresholds.ampl);
%             decision = DMRA1(kf, thresholds);
% %             decision = CheckDecision(decision, kf, kfEthalon, modulations);
%             if (decision == decRight(k))
%                 decRightNum = decRightNum + 1;
%             end
%         end
%         pRight(k, i) = decRightNum / expNum;
%         iteration = iteration + 1;
%         waitbar(iteration / cyclesNum);
%     end
% end
% toc
% close(h);
% 
% figure(2);
% for i = 1 : sigsNum
%     plot(SNR, pRight(i,:), 'marker', mType(i), 'markersize', 10, 'linewidth', 2);
% %     plot(SNR, pRight(i,:), 'linewidth', 2);
%     hold on;
% end
% grid on;
% title('DMRA1'); xlabel('SNR, dB'); ylabel('Probability of right decision');
% legend(decRight, 'location', 'northwest'); legend('show');
% set(gcf, 'color', 'w'); set(groot, 'DefaultAxesFontSize', 18);

%% Offset Test
% modulation = "PSK2";
% idxRight = find(decRight == modulation);
% 
% fOff = linspace(-100, 100, 21);
% % fOff = linspace(-bandsHz(idxRight)/4, bandsHz(idxRight)/4, 41);
% lenOff = length(fOff);
% fOffset = zeros(lenOff, lenSignal);
% tSig = (0 : lenSignal-1) / fs;
% for i = 1 : lenOff
%     fOffset(i, :) = conj(offset) .* exp(1i * 2*pi*(-fOff(i)) * tSig);
% end
% 
% env1 = hilbert(signals(idxRight, :));
% 
% snr = 15;
% expNum = 100;
% cyclesNum = lenOff * expNum;
% 
% prob=zeros(sigsNum, lenOff);
% iteration = 0;
% h = waitbar(0, 'Computing...');
% for i = 1 : lenOff
%     envi = env1 .* fOffset(i,:);
%     for j = 1 : expNum
% %         pos = mod(j, fNum);
%         pos = floor(rand() * (fNum-1));
%         env = envi(pos*lenF+1 : (pos+1)*lenF);
%         env = awgn(env, snr, 'measured');
%         kf = KeyFeatures(env, thresholds.ampl);
%         decision = DMRA1(kf, thresholds);
%         idx = find(decRight == decision);
%         prob(idx, i) = prob(idx, i) + 1;
%         iteration = iteration + 1;
%         waitbar(iteration / cyclesNum);
%     end
% end
% close(h);
% 
% prob = prob / expNum;
% 
% figure(7);
% for i = 1 : sigsNum
%     if i == idxRight
%         line = '-';
%     else
%         line = '--';
%     end
%     plot(fOff, prob(i,:), line, 'marker', mType(i), 'markersize', 10, 'linewidth', 2);
% %     plot(fOff/bandsHz(idxRight), prob(i,:), line, 'marker', mType(i), 'markersize', 10, 'linewidth', 2);
%     hold on;
% end
% grid on;
% title(strcat(modulation, ", ", num2str(snr), " dB"));
% ylabel('Probability of right decision');
% xlabel('\Delta f, Hz');
% % xlabel('$\displaystyle\frac{\delta f}{F}$','interpreter','latex');
% legend(decRight); legend('show');
% set(groot, 'DefaultAxesFontSize', 18); set(gcf, 'color', 'w');




