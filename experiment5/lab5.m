%% Problem 5.1 Impulse Invariance Method and Bilinear Transformation

%%% Butterworth filter

N = 4;
w_cutoff = 0.3 * pi;
T_d = 2;

% impulse invariance
Omega_cutoff = w_cutoff/T_d;

k = 0:N-1;
poly_sk = Omega_cutoff .* exp(1i*pi/(2*N)*(2*k + N + 1));
root_sk = poly(poly_sk);

num = Omega_cutoff^N;
denom = root_sk;


%% Problem 5.2 Characteristics of IIR Filters
close all
clear

f_s = 200;
Rp = 1; % [dB]
Wp = 32/f_s; % [rad/sample]
Ws = 38/f_s; % [rad/sample]
Rs = 25; % [dB]

[N_butter, Wn_butter] = buttord(Wp, Ws, Rp, Rs)
[N_cheb1, Wn_cheb1] = cheb1ord(Wp, Ws, Rp, Rs)
[N_ellip, Wn_ellip] = ellipord(Wp, Ws, Rp, Rs)

[b_butter, a_butter] = butter(N_butter, Wn_butter, 'low'); 
[b_cheb1, a_cheb1] = cheby1(N_cheb1, Rp, Wp);
[b_ellip, a_ellip] = ellip(N_ellip, Rp, Rs, Wp);

[H_butter, w_butter] = freqs(b_butter, a_butter);
[H_cheb1, w_cheb1] = freqs(b_cheb1, a_cheb1);
[H_ellip, w_ellip] = freqs(b_ellip, a_ellip);

%%% Plot the frequency response

figure
set(gcf,'Position',[100 100 800 500])
plot(linspace(0,1,length(H_butter)), abs(H_butter))
title('The frequency response of the analog Butterworth filter')
hold on
grid on
xlabel('n')
ylabel('FFT\{ h(n) \}')
saveas(gcf,'figures/imp_butter.png')

figure
set(gcf,'Position',[100 100 800 500])
plot(linspace(0,1,length(H_cheb1)), abs(H_cheb1));
title('The frequency response of the analog Chebyshev Type I filter')
hold on
grid on
xlabel('n')
ylabel('FFT\{ h(n) \}')
saveas(gcf,'figures/imp_cheb1.png')

figure
set(gcf,'Position',[100 100 800 500])
plot(linspace(0,1,length(H_ellip)), abs(H_ellip))
title('The frequency response of the analog Epilleptic filter')
hold on
grid on
xlabel('n')
ylabel('FFT\{ h(n) \}')
saveas(gcf,'figures/imp_ellip.png')

%%% Plot the frequency response


hn_butter = ifft(H_butter);
hn_cheb1 = ifft(H_cheb1);
hn_ellip = ifft(H_ellip);

figure
subplot(2,1,1)
stem(abs(hn_butter(1:30)))
hold on
grid on
title('Magnitude of the impulse response of the analog Butterworth filter')
xlabel('n')
ylabel('$\abs\{ h(n) \}$', 'interpreter','latex')

subplot(2,1,2)
stem(abs(hn_butter), 'r')
hold on
grid on
title('Imaginary part of the impulse response of the analog Butterworth filter')
xlabel('n')
ylabel('$\mathit{Im}\{ h(n) \}$', 'interpreter','latex')


figure
subplot(2,1,1)
stem(real(hn_cheb1(1:20)))
hold on
grid on
title('Real part of the impulse response of the analog Chebyshev type I filter')
xlabel('n')
ylabel('$\mathit{Re}\{ h(n) \}$', 'interpreter','latex')

subplot(2,1,2)
stem(imag(hn_cheb1(1:20)), 'r')
hold on
grid on
title('Imaginary part of the impulse response of the analog Chebyshev type I filter')
xlabel('n')
ylabel('$\mathit{Im}\{ h(n) \}$', 'interpreter','latex')


figure
subplot(2,1,1)
stem(real(hn_ellip(1:20)))
hold on
grid on
title('Real part of the impulse response of the analog Elliptic filter')
xlabel('n')
ylabel('$\mathit{Re}\{ h(n) \}$', 'interpreter','latex')

subplot(2,1,2)
stem(imag(hn_ellip(1:20)), 'r')
hold on
grid on
title('Imaginary part of the impulse response of the analog Elliptic filter')
xlabel('n')
ylabel('$\mathit{Im}\{ h(n) \}$', 'interpreter','latex')

% Optimum FIR filter

dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1) 10^(-Rs/20)];
[n,fo,ao,w] = firpmord([Wp*f_s Ws*f_s], [1 0], dev, f_s);
b = firpm(n,fo,ao,w);

figure
freqz(b,1)
title('Lowpass Filter Designed to Specifications')




