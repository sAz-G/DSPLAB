function [S] = calc_STFT(x, w, Ts, K)

L = length(w);
N = floor((length(x)-L)/Ts);
S = zeros(N,K);

for i = 1:N
    n_range = (i-1)*Ts+1:(i-1)*Ts+L;
    x_m = x(n_range);
    STFT = fft(w.*x_m, K);

    S(i,:) = STFT;
end

end

