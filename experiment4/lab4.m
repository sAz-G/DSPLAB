%% Problem 7.2 Kaiser Window
clear;
close all; 

w_p = 0.12 * pi;
w_s = 0.18 * pi;
delta_w = (w_s - w_p);

A = 30;

kaiser_window = kwin(delta_w, A).';

N = length(kaiser_window);

fir_lp = firlp((w_s + w_p)/2, N);
kaiser_lp = fir_lp .* kaiser_window;

figure(1)
subplot(2,1,1)
stem(1:N, kaiser_lp)
title("Impulse response of the low-pass")
grid on
hold on
xlabel("n")
ylabel("h(n)")

[H,~] = freqz(kaiser_lp,1);

subplot(2,1,2)

semilogy((1:length(H))/length(H), abs(H))
title("Frequency response of the low-pass")
grid on
hold on
xlabel("\omega/\pi")
ylabel("|H(e^{j\omega})|")
saveas(gcf, "figures/low_pass_kaiser.png")

w_c = 0.15 * pi;

fir_lp = firlp(w_c, N);
kaiser_lp = fir_lp .* kaiser_window;
dirac = zeros(N, 1);
dirac((N+1)/2) = 1;
fft_hp = fft(dirac) - fft(kaiser_lp);
kaiser_hp = ifft(fft_hp);

figure(2)
subplot(2,1,1)
stem(1:N, kaiser_hp)
title("Impulse response of the high-pass")
grid on
hold on
xlabel("n")
ylabel("h(n)")

[H,~] = freqz(kaiser_hp,1);

subplot(2,1,2)

semilogy((1:length(H))/length(H), abs(H))
title("Frequency response of the high-pass")
grid on
xlabel("\omega/\pi")
ylabel("|H(e^{j\omega})|")
saveas(gcf, "figures/high_pass_kaiser.png")

delta_w = 0.25 * pi;
w_0 = 0.3 * pi;

n = -(N-1)/2:(N-1)/2;
h_bp = 2 * (sin(delta_w/2 * n')) ./(pi*n') .* cos(w_0 * n');
h_bp((N+1)/2) = delta_w/pi;
 
kaiser_bp = h_bp .* kaiser_window;

figure(3)
subplot(2,1,1)
stem(1:N, kaiser_bp)
title("Impulse response of the bandpass")
grid on
hold on
xlabel("n")
ylabel("h(n)")

[H,W] = freqz(kaiser_bp,1);

subplot(2,1,2)
semilogy((1:length(H))/length(H), abs(H)) % calculate lim n->0
title("Frequency response of the bandpass")
grid on
xlabel("\omega/\pi")
ylabel("|H(e^{j\omega})|")
saveas(gcf, "figures/band_pass_kaiser.png")

%% Problem 7.5 Hands-on Example: Outdoor Recording

load noise_5insx.mat
load speech.mat

f_s = 8192;

P_signal = sum(abs(x_speech.^2));
P_noise = sum(abs(x_5insx.^2));

SNR_desired = 0;
P_noise_desired = P_signal/(10^(SNR_desired/10));
normalized_noise = (x_5insx/sqrt(P_noise));
new_noise = normalized_noise.*sqrt(P_noise_desired);

signal_1 = x_speech + new_noise;

SNR_1 = 20*log10(norm(x_speech)/norm(new_noise))

SNR_desired = 3;
P_noise_desired = P_signal/(10^(SNR_desired/10));
normalized_noise = (x_5insx/sqrt(P_noise));
new_noise = normalized_noise.*sqrt(P_noise_desired);

signal_2 = x_speech + new_noise;

SNR_2 = 20*log10(norm(x_speech)/norm(new_noise))

N = length(signal_1);

figure
set(gcf,'Position',[100 100 800 400])
subplot(1,2,1)
plot(linspace(-1,1,N),  abs(fftshift(fft(signal_1))))
hold on
grid on
ylim([0, 140])
xlabel("\omega/\pi")
ylabel("|FFT\{x_1(n)\}|")
title("Added signals at 0dB")

subplot(1,2,2)
plot(linspace(-1,1,N), abs(fftshift(fft(signal_2))))
hold on
grid on
ylim([0, 140])
xlabel("\omega/\pi")
ylabel("|FFT\{x_2(n)\}|")
title("Added signals at 3dB")

saveas(gcf, "figures/added_signals.png")

delta_w = 0.1 * pi;
w_cutoff = 0.5 * pi;

A = 30;
kaiser_window = kwin(delta_w, A).';

N_kaiser = length(kaiser_window);

fir_lp = firlp(w_cutoff, N_kaiser);
kaiser_lp = fir_lp .* kaiser_window;


figure
subplot(2,1,1)
stem(1:length(kaiser_lp), kaiser_lp)
title("Impulse response of the low-pass")
grid on
hold on
xlabel("n")
ylabel("h(n)")

[H,~] = freqz(kaiser_lp,1);

subplot(2,1,2)

plot((1:length(H))/length(H), abs(H))
title("Frequency response of the low-pass")
grid on
hold on
xlabel("\omega/\pi")
ylabel("|H(e^{j\omega})|")
saveas(gcf, "figures/low_pass_kaiser.png")

M = length(kaiser_lp) + length(signal_1) - 1;
filterd_signal_1 = ifft(fft(signal_1,M).*fft(kaiser_lp,M));
filterd_signal_2 = ifft(fft(signal_2,M).*fft(kaiser_lp,M));

N = length(filterd_signal_1);

figure
set(gcf,'Position',[100 100 800 400])
subplot(1,2,1)
plot(linspace(-1,1,N),  abs(fftshift(fft(filterd_signal_1))))
hold on
grid on
ylim([0, 140])
xlabel("\omega/\pi")
ylabel("|FFT\{x_1(n)\}|")
title("Added signals at 0dB")

subplot(1,2,2)
plot(linspace(-1,1,N), abs(fftshift(fft(filterd_signal_2))))
hold on
grid on
ylim([0, 140])
xlabel("\omega/\pi")
ylabel("|FFT\{x_2(n)\}|")
title("Added signals at 3dB")

saveas(gcf, "figures/filered_signals.png")

