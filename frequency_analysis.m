% frequency_analysis.m
% BITS F327 AI for Robotics, Hands-On Project 01
%
% Frequency analysis of the teaching signal and the learned signals.
% Reads the recorded data exported from the Serial Monitor.
%
% Columns:
%   1 t          6 alpha1
%   2 theta_d    7 phi1
%   3 theta_hat  8 x
%   4 theta_fix  9 y
%   5 freq      10 freq_hopf

clear; clc; close all;

FILENAME = ['walk_rough.xlsx'];
FMAX     = 10;

D = readmatrix(FILENAME);

t      = D(:,1);
thetaD = D(:,2);
thetaH = D(:,3);
thetaF = D(:,4);
freq   = D(:,5);
alpha1 = D(:,6);
phi1   = D(:,7);
hx     = D(:,8);
hy     = D(:,9);
freqH  = D(:,10);

keep = t >= t(1) + 2;
t = t(keep) - t(find(keep,1));
thetaD = thetaD(keep); thetaH = thetaH(keep); thetaF = thetaF(keep);
freq = freq(keep); alpha1 = alpha1(keep); phi1 = phi1(keep);
hx = hx(keep); hy = hy(keep); freqH = freqH(keep);

N  = length(t);
Fs = (N-1)/(t(end)-t(1));

fprintf('Samples      : %d\n', N);
fprintf('Sample rate  : %.2f Hz\n', Fs);

%% Time domain
figure;
subplot(2,1,1);
plot(t, thetaD, 'k'); hold on;
plot(t, thetaH, 'r', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Amplitude');
legend('Teaching signal', 'Learned signal');
title('Teaching signal and learned signal');
grid on; xlim([0 min(15, t(end))]);

subplot(2,1,2);
plot(t, thetaD, 'k'); hold on;
plot(t, thetaH, 'r', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Amplitude');
title('After convergence');
grid on; xlim([t(end)-10 t(end)]);

%% Spectra
wnd = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
sc  = N/sum(wnd);

Ad = ampSpectrum(thetaD, wnd, sc);
Ah = ampSpectrum(thetaH, wnd, sc);
Af = ampSpectrum(thetaF, wnd, sc);
Ax = ampSpectrum(hx,     wnd, sc);

P = floor(N/2)+1;
f = Fs*(0:P-1)'/N;
band = f <= FMAX;

figure;
plot(f(band), Ad(band), 'k', 'LineWidth', 1.2); hold on;
plot(f(band), Ah(band), 'r', 'LineWidth', 1.2);
plot(f(band), Af(band), 'b--');
plot(f(band), Ax(band), 'g-.');
xlabel('Frequency (Hz)'); ylabel('Amplitude');
legend('Teaching signal', 'Adaptive (harmonic pool)', ...
       'Fixed frequency', 'Hopf oscillator (x)');
title('Single-sided amplitude spectrum');
grid on;

msk = f > 0.3 & f <= FMAX;
[~, id] = max(Ad.*msk);
[~, ih] = max(Ah.*msk);

fprintf('\nDominant frequency, teaching signal   : %.3f Hz\n', f(id));
fprintf('Dominant frequency, learned signal    : %.3f Hz\n', f(ih));
fprintf('Final learned frequency, pool         : %.3f Hz\n', freq(end));
fprintf('Final learned frequency, Hopf         : %.3f Hz\n', freqH(end));
fprintf('Error                                 : %.2f %%\n', ...
        100*abs(f(ih)-f(id))/f(id));

%% Learned frequency
figure;
plot(t, freq, 'r', 'LineWidth', 1.5); hold on;
plot(t, freqH, 'b', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
legend('Harmonic pool', 'Hopf oscillator');
title('Evolution of the learned frequency');
grid on;

%% Tracking error
win = round(Fs);
eH = sqrt(movmean((thetaD-thetaH).^2, win));
eF = sqrt(movmean((thetaD-thetaF).^2, win));

figure;
plot(t, eH, 'r', 'LineWidth', 1.3); hold on;
plot(t, eF, 'b', 'LineWidth', 1.3);
xlabel('Time (s)'); ylabel('RMS error');
legend('Adaptive frequency', 'Fixed frequency');
title('Tracking error comparison');
grid on;

fprintf('\nMean RMS error, adaptive : %.4f\n', mean(eH));
fprintf('Mean RMS error, fixed    : %.4f\n', mean(eF));

%% Phase plane, Cartesian and polar views of the same state
last = t >= t(end)-10;

figure;
subplot(1,2,1);
plot(hx(last), hy(last), 'k'); axis equal; grid on;
xlabel('x'); ylabel('y');
title('Hopf oscillator limit cycle (Cartesian)');

subplot(1,2,2);
rH  = sqrt(hx.^2 + hy.^2);
phH = atan2(hy, hx);
polarplot(phH(last), rH(last), 'k');
title('Same state in polar coordinates');

fprintf('\nMean limit cycle radius : %.3f (expected sqrt(mu) = 1)\n', mean(rH(last)));

%% Pool fundamental in Cartesian form
px = alpha1 .* cos(phi1);
py = alpha1 .* sin(phi1);

figure;
plot(px(last), py(last), 'r'); axis equal; grid on;
xlabel('\alpha_1 cos\phi_1'); ylabel('\alpha_1 sin\phi_1');
title('Fundamental oscillator of the pool, Cartesian view');


function A = ampSpectrum(sig, wnd, sc)
    sig = sig - mean(sig);
    N   = length(sig);
    Y   = fft(sig .* wnd);
    P   = floor(N/2)+1;
    A   = abs(Y(1:P))/N*2*sc;
end
