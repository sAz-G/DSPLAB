%% Problem 6.2 Musicall Tone Systhesis
clear 
close all

f_1 = 392;
f_s = 8192;

N = f_s; % duration = N/f_s for 1s duration
n = 0:(N-1);

y_harmonic = zeros(1,N);
y = sin(2*pi * (f_1/f_s) * n);

for k = 1:8
    p_k = 0.25^(k-1);
    y_harmonic = y_harmonic + p_k*sin(2*pi * (k*f_1/f_s) * n);
end

d = 0.25;
n = [240,7200];
enve = envelope(N,d,n);

tone = y_harmonic.*enve;

figure(1);

hold on
subplot(3,2,1);
plot(1:N, y);
grid on
title("single sinusoidal tone (f_1 = 392 Hz)")
xlabel("n")
ylabel("x(n)")

subplot(3,2,3);
plot(1:N, y_harmonic);
grid on
title("sinusoidal tone with 7 harmonics")
xlabel("n")
ylabel("x(n)")

subplot(3,2,5);
plot((1:N)/N, tone);
grid on
title("sinusoidal tone with 7 harmonics and modified with the envelope")
xlabel("s")
ylabel("x(n)")

subplot(3,2,2);
plot((1:N)-N/2,abs(fftshift(fft(y))));
grid on
title("FFT of the single sinusoidal tone (f_1 = 392 Hz)")
xlabel("Hz")
ylabel('|FFT{x(n)}|')

subplot(3,2,4);
plot((1:N)-N/2, abs(fftshift(fft(y_harmonic))));
grid on
title("FFT of the sinusoidal tone with all harmonics")
xlabel("Hz")
ylabel('|FFT{x(n)}|')

subplot(3,2,6);
plot((1:N)-N/2, abs(fftshift(fft(tone))));
grid on
title("FFt of the modulated sinusodial tone")
xlabel("Hz")
ylabel('|FFT{x(n)}|')

load pianoG3.mat

figure(2)
hold on
subplot(2,1,1);
plot(1:length(g), g);
grid on
title("the recorded tone in the time domain")
xlabel("n")
ylabel("x(n)")

subplot(2,1,2);
plot((1:10162)-10162/2,abs(fftshift(fft(g))));
grid on
title("FFt of the recorded tone")
ylabel("|FFT{x(n)}|")
xlabel("Hz")

y2_harmonic = zeros(1,N);
n = 0:(N-1);

for k = 1:16
    p_k = 0.5^(k-1);
    y2_harmonic = y2_harmonic + p_k*sin(2*pi * (k*f_1/f_s) * n);
end

saveas(gcf, "tone_recorded.png")

d = 0.25;
n = [240,7200];
enve = envelope(N,d,n);

y2_harmonic = y2_harmonic.*enve;

figure(3)
hold on
subplot(2,1,1);
plot((1:length(y2_harmonic))/fs, y2_harmonic);
grid on
title("the recorded tone in the time domain")
xlabel("t [s]")
ylabel("x(n)")


subplot(2,1,2);
plot(1:length(y2_harmonic),abs(fft(y2_harmonic)));
grid on
title("FFt of the recorded tone")
ylabel("|FFT{x(n)}|")
xlabel("Hz")

%% Problem 6.3 Discrete-Time Systems
close all
clear

n = -20:50;
impulse = zeros(1, length(n));
impulse(n == 0) = 1;

b = [0.16 0.48 0.48 0.16];
a = [1 0.13 0.52 0.3];

poles_a = roots(a);
zero_b = roots(b);

abs_poles_a = abs(poles_a)
abs_zeros_b = abs(zero_b)

h = filter(b, a, impulse);
[H,w] = freqz(b,a);

figure(4)
hold on

subplot(2,1,1)
stem(n, h)
title("Impulse response")
ylabel("h(n)")
xlabel("n")
grid on

subplot(2,1,2)
plot(w/pi,20*log10(abs(H)))
title("Frequency response")
ylabel("FFT\{h(n)\}")
xlabel("\omega/\pi")
grid on

saveas(gcf, "response.png")

N = 256;
K = 100;

n2 = 0:N-1;
x_k = zeros(K+1,N);
y_k = zeros(K+1,N);

for k = 0:K
    w_k = pi*k/100;
    x_k(k+1,:) = cos(w_k*n2);
    y_k(k+1, :) = filter(b, a, x_k(k+1,:));
end

y = y_k(:,31:end);
n3 = n2(:,31:end);

figure(5)

hold on
grid on

figure(6)
hold on
grid on

figure(7)
hold on
grid on

for k = 0:10:K
    figure(5)
    plot(n3, abs(y(k+1,:)))
    figure(6)
    plot(n3/256-30/256-0.5*(256 -30)/256, abs(fftshift(fft(y(k+1,:)))))
    figure(7)
    plot(k,mean(abs(y(k+1,:))),'o')
end

figure(5)
title("The magnitude of the filter output")
xlabel("n")
ylabel("|y_k(n)|")
legend("k = 0", "k = 10", "k = 20", "k = 30", "k = 40", "k = 50", "k = 60", "k = 70", "k = 80", "k = 90", "k = 100")
saveas(gcf, "output_time.png")

figure(6)
title("The filter output in the frequency domain")
xlabel("\omega/\pi")
ylabel("|FFT\{y_k(n)\}|")
legend("k = 0", "k = 10", "k = 20", "k = 30", "k = 40", "k = 50", "k = 60", "k = 70", "k = 80", "k = 90", "k = 100")
saveas(gcf, "output_fft.png")

