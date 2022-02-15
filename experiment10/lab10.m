clear
close all

data_14exp48 = read_radar_data('14exp48.dat');
data_14exp49 = read_radar_data('14exp49.dat');
data_105exp71 = read_radar_data('105exp71.dat');
data_105exp94 = read_radar_data('105exp94.dat');

N_14 = length(data_14exp48);
N_105 = length(data_105exp94);

L = 256; % Length of the hamming window
Ts = 10; % Time step of 10 samples
K = 2048; % Frequency bits

w_ham = hamming(L);

stft_14exp48 = calc_STFT(data_14exp48, w_ham, Ts, K).';
stft_14exp49 = calc_STFT(data_14exp49, w_ham, Ts, K).';
stft_105exp71 = calc_STFT(data_105exp71, w_ham, Ts, K).';
stft_105exp94 = calc_STFT(data_105exp94, w_ham, Ts, K).';

y_toward = [0 1];
y_away = [-1 0];

figure

subplot(2,2,1)
samples = [0 length(stft_14exp48(1,:))];
imagesc(samples, y_toward, 10*log(abs(stft_14exp48(:, 500:2500)).^2))
xlabel('samples n')
title('walking toward the camera')
caxis([get_noise_thr(data_14exp48) 150])
c = colorbar;
c.Label.String = 'log-spectrum of the STFT in dB';

subplot(2,2,2)
samples = [0 length(stft_14exp49(1,:))];
imagesc(samples, y_away, 10*log(abs(stft_14exp49(:, 500:2500)).^2))
xlabel('samples n')
title('walking away the camera')
caxis([get_noise_thr(data_14exp49) 150])
c = colorbar;
c.Label.String = 'log-spectrum of the STFT in dB';

subplot(2,2,3)
samples = [0 length(stft_105exp71(1,:))];
imagesc(samples, y_away, 10*log(abs(stft_105exp71(:, 500:2500)).^2))
xlabel('samples n')
ylabel('frequency \omega/\pi')
title('walking away the camera')
caxis([get_noise_thr(data_105exp71) 150])
c = colorbar;
c.Label.String = 'Log-Spectrum of the STFT in dB';

subplot(2,2,4)
samples = [0 length(stft_105exp94(1,:))];
imagesc(samples, y_toward, 10*log(abs(stft_105exp94(:, 500:2500)).^2))
xlabel('samples n')
ylabel('frequency \omega/\pi')
title('walking toward the camera')
caxis([get_noise_thr(data_105exp94) 150])
c = colorbar;
c.Label.String = 'log-spectrum of the STFT in dB';
