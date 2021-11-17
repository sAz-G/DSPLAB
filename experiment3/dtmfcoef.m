function [fr, fc] = dtmfcoef(dtft_dtmf_signal,fs)
%DTMFCOEF Computes the row and collum frequency in Hz
%   [fr, fc] = dtmfcoef(dtft_dtmf_signal,fs)
%   
%   dtft_dtmf_signal = the approximated DTFT of the DTMF signal
%   fs = sample rate

    absolute_signal = abs(dtft_dtmf_signal); % Absolute value of the DTMF signal

    local_max = islocalmax(absolute_signal); % Safe all local maxima
    position = 1:length(absolute_signal); % Position vector
    position_of_max = position(local_max);  % Position of the local maxima
    value_of_max = absolute_signal(local_max); % Value of the local maxima

    % Extract the first two (of four) main peaks
    [~, pos_four_largest] = maxk(abs(value_of_max), 4); % Search for the 4 largest maxima
    pos_first_two = mink(pos_four_largest, 2); %  Get the position of the first two main peaks

    % Get the position and the corresponding frequencies
    maxima_frequency = position_of_max(pos_first_two); % Safe the position of the first two maxima
    normalized_frequency = maxima_frequency/length(dtft_dtmf_signal); % Determine the normalized angular frequency of the peaks
    symbol_freq = normalized_frequency*fs; % Calculate the frequency in Hz
    fr = min(symbol_freq); 
    fc = max(symbol_freq);
end

