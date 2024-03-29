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

w = linspace(0,1,1000);

[N_butter, Wn_butter] = buttord(Wp, Ws, Rp, Rs, 's')
[N_cheb1, Wn_cheb1] = cheb1ord(Wp, Ws, Rp, Rs, 's')
[N_ellip, Wn_ellip] = ellipord(Wp, Ws, Rp, Rs, 's')

[b_butter, a_butter] = butter(N_butter, Wn_butter,'s'); 
[b_cheb1, a_cheb1] = cheby1(N_cheb1, Rp, Wp,'s');
[b_ellip, a_ellip] = ellip(N_ellip, Rp, Rs, Wp,'s');

sys_butter = tf(b_butter, a_butter);
sys_cheb1 = tf(b_cheb1, a_cheb1);
sys_ellip = tf(b_ellip, a_ellip);

[H_butter, w_butter] = freqs(b_butter, a_butter, w);
[H_cheb1, w_cheb1] = freqs(b_cheb1, a_cheb1, w);
[H_ellip, w_ellip] = freqs(b_ellip, a_ellip, w);

%%% Plot the frequency response
figure
set(gcf,'Position',[100 100 800 700])
subplot(3,1,1)
plot(w, mag2db(abs(H_butter)))
title('The frequency response of the analog Butterworth filter')
hold on
grid on
xlabel('f/f_s')
ylabel('FFT\{h(t)\} [dB]')

subplot(3,1,2)
plot(w, mag2db(abs(H_cheb1)));
title('The frequency response of the analog Chebyshev Type I filter')
hold on
grid on
xlabel('f/f_s')
ylabel('FFT\{h(t)\} [dB]')

subplot(3,1,3)
plot(w*f_s, mag2db(abs(H_ellip)))
title('The frequency response of the analog Elliptic filter')
hold on
grid on
xlabel('f/f_s')
ylabel('FFT\{h(t)\} [dB]')

saveas(gcf,'figures/analog_frequency_response.png')

%%% Plot the impulse response

[imp_butter, t_butter] = impulse(sys_butter);
[imp_cheb1, t_cheb1] = impulse(sys_cheb1);
[imp_ellip, t_ellip] = impulse(sys_ellip);

figure
set(gcf,'Position',[100 100 800 700])
subplot(3,1,1)
plot(t_butter, imp_butter)
hold on
grid on
title('The impulse response of the analog Butterworth filter')
xlabel('t [s]')
ylabel('h(t)')
xlim([0 400])
ylim([-0.04 0.06])

subplot(3,1,2)
plot(t_cheb1, imp_cheb1)
hold on
grid on
title('The impulse response of the analog Chebyshev type I filter')
xlabel('t [s]')
ylabel('h(t)')
xlim([0 400])
ylim([-0.04 0.06])

subplot(3,1,3)
plot(t_ellip, imp_ellip)
hold on
grid on
title('The impulse response of the analog Elliptic filter')
xlabel('t [s]')
ylabel('h(t)')
xlim([0 400])
ylim([-0.04 0.06])

saveas(gcf,'figures/analog_impulse_response.png')

%%% zero ploe plot

figure
zplane(b_butter, a_butter)
title('Pole-Zero Plot of the Butterworth filter')
grid on
saveas(gcf,'figures/butter_zero_pole.png')

figure
zplane(b_cheb1, a_cheb1)
title('Pole-Zero plot of the Chebyshev type I filter')
grid on
saveas(gcf,'figures/cheb1_zero_pole.png')

figure
zplane(b_ellip, a_ellip)
title('Pole-Zero plot of the Elliptic filter')
grid on
saveas(gcf,'figures/ellip_zero_pole.png')


%%% Optimum FIR filter

dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1) 10^(-Rs/20)];
[n,fo,ao,w] = firpmord([Wp*f_s Ws*f_s], [1 0], dev, f_s);
b = firpm(n,fo,ao,w);
[H, w] = freqz(b,1,1024,f_s);

figure
plot(w, mag2db(abs(H)))
title('Lowpass Filter Designed to Specifications')
grid on
xlabel('f [Hz]')
ylabel('FFT\{h(n)\} [dB]')
saveas(gcf,'figures/optimum_FIR.png')

%% Problem 5.3: IIR Filtering of Sinusoids

clear
close all

%%% Part 1

% Design the filter

clear
close all

Wp_band = [125 175]; % [Hz]
Ws_band = [115 185]; % [Hz]
Rp_band = 1; % [dB]
Rs_band = 40; % [dB]

[N_cheb1, Wn_cheb1] = cheb1ord(Wp_band, Ws_band, Rp_band, Rs_band, 's');
f_s_band = 1000; % [Hz]

[z,p,k] = cheby1(N_cheb1,1,2*pi*Wn_cheb1,'bandpass','s');

[num,den] = zp2tf(z,p,k);
[h,w] = freqs(num,den);

figure
plot(w/(2*pi), mag2db(abs(h)))
hold on
xlim([0 500])
grid
legend('Magnitude response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')

[numd,dend] = bilinear(num,den,f_s_band,150);

fvtool(numd,dend,'Fs',f_s_band)

figure
subplot(2,1,1)
plot(w/(2*pi), mag2db(abs(h)))
title('The frequency response of the analog bandpass filter')
xlim([1 500])
hold on
grid on
xlabel('f [Hz]')
ylabel('FFT\{h(t)\} [dB]')

[hd, wd] = freqz(numd, dend);

subplot(2,1,2)
plot(wd, mag2db(abs(hd)))
title('The frequency response of the digital bandpass filter')
hold on
grid on
xlabel('\omega [rad/sample]')
ylabel('FFT\{h(n)\} [dB]')

saveas(gcf,'figures/FIR_bandpasss.png')


%%% Filter the sequence
T_s_band = 1/f_s_band;
samples = 300;
t = linspace(0,samples*T_s_band, samples+1);
x_t = 5*sin(200*pi*t)+2*cos(300*pi*t);

y_t = filter(numd, dend, x_t);

L = 1000;
fft_xt = fft(x_t, L);
% fft_xt = fft_xt(1:L/2);

fft_yt = fft(y_t, L);
% fft_yt = fft_yt(1:L/2);

f = linspace(-f_s_band,f_s_band,L);

figure
subplot(2,1,1)
plot(f, fftshift(abs(fft_xt)))
xlabel("f [Hz]")
ylabel('FFT\{x(t)\}')
title("Fourier transform of the signal")
legend("FFT of x(t)")
grid on

subplot(2,1,2)
plot(f, fftshift(abs(fft_yt)))
xlabel("f [Hz]")
ylabel('FFT\{y(t)\}')
title("Fourier transform of the filtered signal")
legend("FFT of y(t)")
grid on

saveas(gcf,'figures/FFT_signal_bandpass.png')

%% Part 2
close all

N_e = 10;

w_c_e = 0.44 * pi;
w_bandwidth_e = 0.08 * pi;

Rp_e = 1;
Rs_e = 100;
Wp_e = [w_c_e - w_bandwidth_e/2, w_c_e + w_bandwidth_e/2];

[b_e, a_e] = ellip(N_e, Rp_e, Rs_e, Wp_e/pi, 'stop');
[h_e, w_e] = freqz(b_e, a_e, 1024);

figure
plot(w_e/pi, mag2db(abs(h_e)))
title('The frequency response of the digital elliptic bandsop filter')
hold on
grid on
xlabel('\omega/\pi [rad/sample]')
ylabel('FFT\{h(n)\} [dB]')

saveas(gcf,'figures/FIR_bandstop.png')

n = 0:200;
x_n_e = sin(0.44 * pi * n);
y_n_e = filter(b_e, a_e, x_n_e);

L = 1000;
f = linspace(-1,1,L);

figure
subplot(2,1,1)
plot(f, fftshift(abs(fft(x_n_e, L))))
xlabel('\omega/\pi [rad/sample]')
ylabel('FFT\{x(n)\}')
title("Fourier transform of the signal")
legend("FFT of x(n)")
grid on

subplot(2,1,2)
plot(f, fftshift(abs(fft(y_n_e, L))))
xlabel('\omega/\pi [rad/sample]')
ylabel('FFT\{y(n)\}')
title("Fourier transform of the filtered signal")
legend("FFT of y(n)")
grid on

saveas(gcf,'figures/FFT_signal_bandstop.png')

figure
subplot(2,1,1)
plot(x_n_e)
xlabel("n")
ylabel('x(n)')
title("Sample sequence in the time-domain")
legend("x(n)")
ylim([-1 1])
grid on

subplot(2,1,2)
plot(y_n_e)
xlabel("n")
ylabel('y(n)')
title("Filtered signal in the time-domain")
ylim([-1 1])
legend("y(n)")
grid on

saveas(gcf,'figures/signal_bandstop.png')

