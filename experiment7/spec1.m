function out_IxxM = spec1(x_n, window, L, resolution)
%SPEC1 Summary of this function goes here
%   x_n = the signal data
%   w the window
%   L the number of segments

out_IxxM = zeros(1, resolution);
M = length(window);

for l = 1:L
    nl = (l-1)*M+1:l*M;
    xl_n = x_n(nl);
    
    noralizing_factor = sum(abs(window).^2);
    
    windowed_xl_n = window .* xl_n;
    windowed_Xl  = fft(windowed_xl_n, resolution);
    IxxM_l =  windowed_Xl .* conj(windowed_Xl)/noralizing_factor;
    out_IxxM = out_IxxM + IxxM_l;
end

out_IxxM = out_IxxM/L;

end
