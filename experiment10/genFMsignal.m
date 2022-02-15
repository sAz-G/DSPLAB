function [signal] = genFMsignal(alpha, beta, gamma, fc, fs, tStart, tStop)

t = tStart:1/fs:tStop;
% wi = 2*pi*fc + 2*pi*alpha*cos(2*pi*beta*t + gamma);
phi = 2*pi*fc*t + alpha/beta*cos(2*pi*beta*t + gamma);

signal = A*cos(phi);

end

