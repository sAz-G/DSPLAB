clear
close all

%% Problem 4.1 Windowed, Averaged Periodograms
load ar7.mat

resolution = 512;

L = [1, 2, 5, 10];
M = length(X)./L;

w=pi*(0:1/resolution*2:1-1/resolution);

for i = 1:4
    window = ones(M(i), 1);
    spectrum = spec1(X,window,L(i), resolution);

    subplot(length(L),2,i*2-1)
    plot(w/pi, 10*log10(spectrum(1:resolution/2)))
    title("Rect: L = " + num2str(L(i)))

    window = hamming(M(i));
    spectrum = spec1(X, window, L(i), resolution);

    subplot(length(L), 2, i*2)
    plot(w/pi, 10*log10(spectrum(1:resolution/2)))
    title("Hamming: L = " + num2str(L(i)))
end

% Log-Spectrum

L2 = 5;
alpha = 0.05;
M2 = length(X)/L2;
window2 = hamming(M2);
LogSpectrum = 10*log10(real(spec1(X, window2, L2, resolution)));

confidence_interval_lb = LogSpectrum - 10*log10(exp(1)) * norminv(0.025)/sqrt(L2);
confidence_interval_ub = LogSpectrum + 10*log10(exp(1)) * norminv(0.025)/sqrt(L2);

figure
plot(w/pi, LogSpectrum(1:resolution/2))
hold on

title("The log-spectrum for L=5 and a hamming window")
plot(w/pi, confidence_interval_lb(1:resolution/2), "r")
plot(w/pi, confidence_interval_ub(1:resolution/2), "r")

% simon.tien@stud.tu-darmstadt.de

%% Problem 4.4 Biomedical Data
clear
close all

dt = 0.01;
f_s = 1/dt;
f_s_filter = 200;

% Filter
N = 30;
f_c = 40;
[b, a] = butter(N, f_c/(f_s_filter/2));

% preprocess pulse signal
[y_pulse, t_pulse] = readfile(['testbloodpulse.txt']);
y_pulse_det = detrend(y_pulse);
y_pulse_filtered = filter(b,a,y_pulse_det);

% preprocess fast respiration signal
[y_fast, t_fast] = readfile('testfast.txt');
y_fast_det = detrend(y_fast);
y_fast_filtered = filter(b,a,y_fast_det);

% preprocess slow respiration signal
[y_slow, t_slow] = readfile('testslow.txt');
y_slow_det = detrend(y_slow);
y_slow_filtered = filter(b,a,y_slow_det);

% Parameters for ploting
offset = 102;

duration_pulse = 10;
duration_resp = 20;

index_plot_pulse = offset:duration_pulse*f_s+offset+1;
index_plot_resp = offset:duration_resp*f_s+offset+3;

% Plot pulse respiration signal
y_pulse_filtered_part = y_pulse_filtered(index_plot_pulse);

figure
plot(t_pulse(index_plot_pulse), y_pulse_filtered_part)
title('Excerpt of 10 s from the pulse measurements')
grid on
xline(5, 'r')
xline(6, 'r')
xlim([t_pulse(offset), t_pulse(duration_pulse*f_s+offset)])
xlabel("s")
ylabel("mV")
legend("pulse","one cardiac cycle")


% Plot slow respiration signal
y_slow_filtered_part = y_slow_filtered(index_plot_resp);

figure
subplot(2,1,1)

plot(t_slow(index_plot_resp), y_slow_filtered_part)
title('Excerpt of 10 s from the slow respiration measurements')
grid on
xline(4.5, 'r')
xline(14, 'r')
xlim([t_slow(offset), t_slow(duration_resp*f_s+offset)])
xlabel("s")
ylabel("mV")
legend("pulse","one cardiac cycle")

% Plot fast respiration signal
y_fast_filtered_part = y_fast_filtered(index_plot_resp);

subplot(2,1,2)
plot(t_fast(index_plot_resp), y_fast_filtered_part)
title('Excerpt of 10 s from the fast respiration measurements')
grid on
xline(6, 'r')
xline(8, 'r')
xlim([t_fast(offset), t_fast(duration_resp*f_s+offset)])
xlabel("s")
ylabel("mV")
legend("pulse","one cardiac cycle")

% Calculate and plot the Periodogram 

L = 1;
spectrum_pulse = spec1(y_pulse_filtered_part.', ones(1, length(y_pulse_filtered_part)/L), L, 512);
spectrum_slow = spec1(y_slow_filtered_part.', ones(1, length(y_slow_filtered_part)/L), L, 512);
spectrum_fast = spec1(y_fast_filtered_part.', ones(1, length(y_fast_filtered_part)/L), L, 512);

L=6;
spectrum_pulse_avg = spec1(y_pulse_filtered_part.', ones(1, length(y_pulse_filtered_part)/L), L, 512);
spectrum_slow_avg = spec1(y_slow_filtered_part.', ones(1, length(y_slow_filtered_part)/L), L, 512);
spectrum_fast_avg = spec1(y_fast_filtered_part.', ones(1, length(y_fast_filtered_part)/L), L, 512);

% Plot pulse respiration signal
figure
semilogy(spectrum_pulse_avg)
title('Spectrum of the pulse measurements')
grid on
xlabel("s")
ylabel("mV")
legend("pulse","one cardiac cycle")


% Plot slow respiration signal
figure
subplot(2,1,1)

semilogy(abs(spectrum_slow_avg))
title('Spectrum of the slow respiration measurements')
grid on
xlabel("1/Hz")
ylabel("mV")
legend("spectrum","one cardiac cycle")

% Plot fast respiration signal
subplot(2,1,2)
semilogy(abs(spectrum_fast_avg))
title('Spectrum of the fast respiration measurements')
grid on
xlabel("1/Hz")
ylabel("mV")
legend("pulse","one cardiac cycle")



