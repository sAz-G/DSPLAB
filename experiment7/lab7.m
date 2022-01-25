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

