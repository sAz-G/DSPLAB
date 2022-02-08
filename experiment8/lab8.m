%% Problem 6.1 AR process

clear
close all

N = 2^10;
X_n = randn(N, 1);
a = [1, 0.5, 0.7, 0.2]; % X(n) + 0.5X(n − 1) + 0.7X(n − 2) + 0.2X(n − 3) = Z(n)

y_n = filter(1,a,X_n); % Filter the static white noise
arcoeffs = aryule(y_n,3); % Estimate the AR coefficients by Yule-Walker equations

[pxx_aryule, w_aryule] = freqz(1,arcoeffs,N);
[pxx_pyulear, w_pyule] = pyulear(y_n,3,N);
[cxx, w] = freqz(1,a,N);

figure
plot(w_aryule/pi, abs(pxx_aryule), "r")
title("Power Spectral Density Esitmate")

hold on
plot(w_pyule/pi, abs(pxx_pyulear), "b")

plot(w/pi, abs(cxx), "k")

grid on
% ylim([-15 10])
ylabel('$|P_{xx}|$','Interpreter','latex')
xlabel('$\omega/pi$','Interpreter','latex')
legend("aryule PSD Estimate", "pyulear PSD Estimate", "true PSD")

% [max_spec, pos_spec] = max(abs(pxx_aryule));
% xline(w(pos_spec)/pi,'r');
% [max_spec, pos_spec] = max(abs(pxx_aryule));
% xline(w(pos_spec)/pi,'b');
% [max_spec, pos_spec] = max(abs(cxx));
% xline(w(pos_spec)/pi,'k');

saveas(gcf,'figures/spectrum_AR.png')

%% Problem 6.3 Order Selection

load arunknown.mat
N = length(X);
M = 10;

aic_m = zeros(10,1);
mdl_m = zeros(10,1);
var_m = zeros(10,1);

for m = 1:M
    a_m = aryule(X,m);
    
%     Y = filter(a_m ,1 ,X);
%     Y = Y(m:N);
    
    Y = zeros(N-m, 1);
    for n = m+1:N
        Y(n) = X(n) + sum(a_m(2:m+1) .* fliplr(X(n-m:n-1)));
    end

    var_Zm = 1/(N - m) * sum(Y.^2);
    
    % Akaike's Information Criterion (AIC)
    aic_m(m) = log(var_Zm) + m * 2/N;

    % Minimum Description length (MDL)
    mdl_m(m) = log(var_Zm) + m * log(N)/N;
    
    var_m(m) = var_Zm;
end

[min_aic, pos_aic] = min(aic_m);
[min_mdl, pos_mdl] = min(mdl_m);

figure
plot(1:M, var_m)
title("The variance $\sigma_{Zm}$",'Interpreter','latex')
grid on
xlabel('m')
ylabel("$\sigma_{Zm}$",Interpreter="latex",FontSize=12,FontWeight="bold")

figure
plot(1:M, aic_m, 'k')
hold on
plot(1:M, mdl_m, 'b')
grid on
plot(pos_aic,aic_m(pos_aic), 'k*')
plot(pos_mdl,mdl_m(pos_mdl), 'b*')

title("The AIC and the MDL over the order m")
xlabel('m')
ylabel('AIC(m) and MDL(m)')
legend("Akaike's Information Criterion (AIC)", 'Minimum Description length (MDL)')

saveas(gcf,'figures/AIC_MDL.png')

%% Problem 6.4 Signal Analysis
clear
close all

load s5.mat
f_s = 8e3;

SH = s5(15600:16300);
AA = s5(16800:17500);

% soundsc(SH, f_s);
% soundsc(AA, f_s);

figure(5)
subplot(2,1,1)
ylim([-2e4 2e4])
plot((15600:16300)/f_s, SH)
title('Signal of the SH segment')
ylabel('signal s5')
xlabel('t [s]')

subplot(2,1,2)
plot((16800:17500)/f_s, AA)
title('Signal of the AA segment')
ylim([-2e4 2e4])
ylabel('signal s5')
xlabel('t [s]')

M = 10;
N = 2^10;

a_SH = aryule(SH, M);
a_AA = aryule(AA, M);

% Order Selection:
N_AA = length(AA);
N_SH = length(SH);
M = 35;

aic_m_aa = zeros(M,1);
aic_m_sh = zeros(M,1);
var_m_aa = zeros(M,1);
var_m_sh = zeros(M,1);

for m = 1:M

    a_m_aa = aryule(AA,m);
    a_m_sh = aryule(SH,m);
    
%     Y_aa = filter(a_m_aa, 1, AA);
%     Y_aa = Y_aa(m:N_AA);
%     
%     Y_sh = filter(a_m_sh, 1, SH);
%     Y_sh = Y_sh(m:N_SH);
    
    Y_aa = zeros(N_AA-m, 1);
    for n = m+1:N_AA
        Y_aa(n) = AA(n) + sum(a_m_aa(2:m+1) .* fliplr(AA(n-m:n-1).'));
    end

    Y_sh = zeros(N_SH-m, 1);
    for n = m+1:N_SH
        Y_sh(n) = SH(n) + sum(a_m_sh(2:m+1) .* fliplr(SH(n-m:n-1).'));
    end

    var_Zm_AA = 1/(N_AA - m) * sum(Y_aa.^2);
    var_Zm_SH = 1/(N_SH - m) * sum(Y_sh.^2);
    
    % Akaike's Information Criterion (AIC)
    aic_m_aa(m) = log(var_Zm_AA) + m * 2/N_AA;
    aic_m_sh(m) = log(var_Zm_SH) + m * 2/N_SH;
    
    var_m_aa(m) = var_Zm_AA;
    var_m_sh(m) = var_Zm_SH;
end

[~, pos_aa] = min(aic_m_aa)
[~, pos_sh] = min(aic_m_sh)

[pxx_SH, ~] = freqz(1, a_SH, N);
pxx_SH = abs(pxx_SH)/max(abs(pxx_SH));
[pxx_AA, w] = freqz(1, a_AA, N);
pxx_AA = abs(pxx_AA)/max(abs(pxx_AA));

L = 8;
[pxx_SH_welch, ~] = pwelch(SH, hamming(N/L), (N/L)/2, N);
pxx_SH_welch = abs(pxx_SH_welch)/max(abs(pxx_SH_welch));
[pxx_AA_welch, ~] = pwelch(AA, hamming(N/L), (N/L)/2, N);
pxx_AA_welch = abs(pxx_AA_welch)/max(abs(pxx_AA_welch));

figure(6)
subplot(2,1,1)
plot(w/pi, abs(pxx_SH), 'k');
hold on
plot(linspace(0,1,length(pxx_SH_welch)), abs(pxx_SH_welch), 'b');
grid on
title('AR(10) spectrum of the SH segment')
xlabel('$\omega/pi$','Interpreter','latex')
ylabel('$|\hat{C}_{xx}|$','Interpreter','latex')
legend('pwelch - PSD Estimator', 'aryule - PSD Estimator')

subplot(2,1,2)
plot(w/pi, abs(pxx_AA), 'k');
hold on
plot(linspace(0,1,length(pxx_AA_welch)), abs(pxx_AA_welch), 'b');
grid on
title('AR(10) spectrum of the AA segment')
xlabel('$\omega/pi$','Interpreter','latex')
ylabel('$|\hat{C}_{xx}|$','Interpreter','latex')
legend('pwelch - PSD Estimator', 'aryule - PSD Estimator',"Location","Northwest")




