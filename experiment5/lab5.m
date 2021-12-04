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
xlabel('\omega/\pi [rad/sample]')
ylabel('FFT\{h(n)\} [dB]')

subplot(3,1,2)
plot(w, mag2db(abs(H_cheb1)));
title('The frequency response of the analog Chebyshev Type I filter')
hold on
grid on
xlabel('\omega/\pi [rad/sample]')
ylabel('FFT\{h(n)\} [dB]')

subplot(3,1,3)
plot(w, mag2db(abs(H_ellip)))
title('The frequency response of the analog Epilleptic filter')
hold on
grid on
xlabel('\omega/\pi [rad/sample]')
ylabel('FFT\{h(n)\} [dB]')

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

subplot(3,1,2)
plot(t_cheb1, imp_cheb1)
hold on
grid on
title('The impulse response of the analog Chebyshev type I filter')
xlabel('t [s]')
ylabel('h(t)')

subplot(3,1,3)
plot(t_ellip, imp_ellip)
hold on
grid on
title('The impulse response of the analog Elliptic filter')
xlabel('t [s]')
ylabel('h(t)')

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
title('Pole-Zero plot of the Epilleptic filter')
grid on
saveas(gcf,'figures/ellip_zero_pole.png')


%%% Optimum FIR filter

dev = [(10^(Rp/20)-1)/(10^(Rp/20)+1) 10^(-Rs/20)];
[n,fo,ao,w] = firpmord([Wp*f_s Ws*f_s], [1 0], dev, f_s);
b = firpm(n,fo,ao,w);

figure
freqz(b,1)
title('Lowpass Filter Designed to Specifications')
saveas(gcf,'figures/optimum_FIR.png')

%% Problem 5.3: IIR Filtering of Sinusoids

clear
close all

%%% Part 1

% Design the filter

f_s_band = 1000; % [Hz]
Wp_band = [125 175]/f_s_band; % [rad/sample]
Ws_band = [115 185]/f_s_band; % [rad/sample]
Rp_band = 1; % [dB]
Rs_band = 40; % [dB]

w = linspace(0,1,1000);

[N_band, Wn_band] = buttord(Wp_band, Ws_band, Rp_band, Rs_band, 's');
[z_band, p_band, k_band] = butter(N_band, Wn_band,'s');

[b_band, a_band] = zp2tf(z_band, p_band, k_band);
H_band = freqs(b_band, a_band, w);

[bd_band, ad_band] = impinvar(b_band, a_band, f_s_band);

Hd_band = freqz(bd_band, ad_band, w);

figure
subplot(2,1,1)
plot(w, mag2db(abs(H_band)))
title('The frequency response of the bandpass filter')
hold on
grid on
xlabel('f [Hz]')
ylabel('FFT\{h(n)\} [dB]')

subplot(2,1,2)
plot(w, mag2db(abs(Hd_band)))
title('The frequency response of the bandpass filter')
hold on
grid on
xlabel('f [Hz]')
ylabel('FFT\{h(n)\} [dB]')


%%% Filter the sequence
T_s_band = 1/f_s_band;
samples = 300;
t = linspace(0,samples*T_s_band, samples+1);
x_t = 5*sin(200*pi*t)+2*cos(300*pi*t);

y_t = filter(b_band, a_band, x_t);

figure
subplot(2,1,1)
plot(t, x_t)

subplot(2,1,2)
plot(t, y_t)



