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
    p_k = 0.4^(k-1);
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

subplot(3,2,3);
plot(1:N, y_harmonic);

subplot(3,2,5);
plot(1:N, tone);

subplot(3,2,2);
plot(1:N,abs(fftshift(fft(y))));

subplot(3,2,4);
plot(1:N,abs(fftshift(fft(y_harmonic))));

subplot(3,2,6);
plot(1:N, abs(fftshift(fft(tone))));

load pianoG3.mat

figure(2)
hold on
subplot(2,1,1);
plot(1:10162, g);

subplot(2,1,2);
plot(1:10162,abs(fftshift(fft(g))));

%% Problem 6.3 Discrete-Time Systems
close all
clear

n = -20:50;
impulse = zeros(1, length(n));
impulse(n == 0) = 1;

a = [0.16 0.48 0.48 0.16];
b = [1 0.13 0.52 0.3];

poles_a = roots(a);
zero_b = roots(b);

abs_poles_a = abs(poles_a)
abs_zeros_b = abs(zero_b)

h = filter(b, a, impulse);
[H,w] = freqz(b,a);

figure(3)
hold on

subplot(2,1,1)
stem(n, h)
grid on

subplot(2,1,2)
plot(w/pi,20*log10(abs(H)))
grid on

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

figure(4)
hold on
grid on
for k = 0:K
    plot(n3, 20*log10(abs(y(k+1,:))))
end

