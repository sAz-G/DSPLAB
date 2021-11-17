function dft_x = dft(x_n)
%DTF Computed the DTF of a given signal x
%   dft_x = dft(x_n)
%
%   x_n = the signal in the time domain

    N = length(x_n); % Length of the signal
    
    % Computation of the DFT
    n =  0:(N-1);
    exp_mat = n'*n;
    exp_fft = exp(-1i*2*pi/N*exp_mat);
    dft_x = exp_fft * x_n;
end

