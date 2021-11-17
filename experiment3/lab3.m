%% Computationat Complexity of DFT/FFT
close all
clear

k = 4:8;
N_row = 2.^k;

time_dft = zeros(length(k), 1);
time_fft = zeros(length(k), 1);

L = 1000;


for i = 1:length(k)
    
    for l = 1:L
        signal = randn(1,max(N_row));
        
        N = N_row(i);
        x = signal(1:N-1);

        tic
        dft_x = dft(x.');
        time_dft(i) = time_dft(i) + toc;

        tic
        fft_x = fft(x.');
        time_fft(i) = time_fft(i) + toc;
    end
    
    time_dft(i) = time_dft(i)/L;
    time_fft(i) = time_fft(i)/L;
end

figure(1)

semilogy(N_row, time_dft)
hold on
grid on
semilogy(N_row, time_fft)

title("runtime of the dft and fft algorithm ")
legend("dft", "fft")
xlabel("length of the signal N")
ylabel("average runtime [s]")

saveas(gcf, "figures/runtime.png")

%% DTMF Signal Decoding
close all

fs = 8000;

dtmf_signal_seven = dtmfdial('7', fs);
dtft_dtmf_signal = fft(dtmf_signal_seven);

figure(2)
plot((1:length(dtft_dtmf_signal))/length(dtft_dtmf_signal)*fs, abs(dtft_dtmf_signal))
hold on
grid on

fc_desired = 1209;
fr_desired = 852;

xline(fr_desired,'r')
xline(fc_desired,'k')
xline(fs-fr_desired,'r')
xline(fs-fc_desired,'k')

title("DFT of the DTMF signal of symbol '7'")
xlabel("\omega/\pi")
ylabel("FFT\{x(n)\}")
legend("DTMF", "f_r", "f_c")

% Estimate frequency of the approximated DTFT of the DTMF signal
[fr_est, fc_est] = dtmfcoef(dtft_dtmf_signal, fs);
dial = get_dtmf_dial(fr_est, fc_est)

load mynumber.mat;

dial_string = repmat(['-'],1,11);
symbol_count = 11;
dtmf_length = 256;

for i = 1:symbol_count
    offset = (i-1)*length(xx)/symbol_count;
    
    range = offset+1:offset+dtmf_length;
    symbol_timedomain = xx(range);
    symbol_fft = fft(symbol_timedomain);
    
    [fr_est, fc_est] = dtmfcoef(symbol_fft, fs);
    dial_est = get_dtmf_dial(fr_est, fc_est);
    dial_string(i) = dial_est;
end

dial_string

