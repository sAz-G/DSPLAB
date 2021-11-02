%% solution for experiment 1

%% Problem_4.1

M = magic(5);

add_M       = sum(M);
add_trans_M = sum(M');

row_1 = M(1,:);
col_3 = M(:,3);

col_1_to_3_of_2_end = M(2:end, 1:3);

greater_10 = find(M>10);
less_4     = find(M<4);

%% Problem 4.2

g = fibonacci(12);
f_ratio = g(2:12)./g(1:11);

figure(1)
grid on
plot(f_ratio)

%% Problem 4.3
clear
x = rand(100000,1);

min_x = min(x);
max_x = max(x);

mn_x = mean(x);
deviation_x = std(x);
var_x       = var(x);

y = 4*x-2;

mean_y      = mean(y);
deviation_y = std(y);

%%%%% do tho randn part 

%% Part 4.4
clc, clear
V_ml = 330;
r_cm = 0.5:.01:10;
A = 2*pi*r_cm.^2+2*pi*(V_ml./(pi*r_cm));
h = V_ml./(pi*r_cm.^2);
opt_r = r_cm(A == min(A));
h_opt = h(r_cm == opt_r);
min_A = min(A);

figure(2)

plot(r_cm, A, 'r')

grid on
hold on 

plot(opt_r, min_A, 'b*')
xlabel('r')
ylabel('A')
title('Dimensions of Common Beverage Can')


%% Part 4.6
clear, clc
n = 0:100;

F_Hz = 1;
T_s = 0.05;

s = sin(2*pi*F_Hz*n*T_s);

S = fft(s,128);
P = S.*conj(S);
w = (0:127)./128;

close all
figure(3)
plot(n, s)
grid on
hold on
stem(n,s)
hold on
plot(n*T_s, s)
legend('discrete', 'continues', 'stem')

figure(4)
plot(2*w, P)
grid on 
hold on 
plot(w/T_s, P)

%%%% disturbance 
s2 = s + sin(2*pi*4*n*T_s);
S2 = fft(s2, 128);

b = [1 1 1 1]/4;
a = 1;

[H, w1] = freqz(b,a);
figure(6)
plot(w1/(2*pi*T_s), abs(H))
grid on 

figure(5)
plot(n, s)
hold on 
plot(n,s2)

%%% filter the signal 
sf_2 = filter(b,a,s2);
SF_2 = fft(sf_2,128);

sf = filter(b,a, s);
SF = fft(sf, 128);

figure(7)
plot(w/T_s, abs(SF_2));
grid on 
hold on
plot(w/T_s, abs(SF));

figure(8)
plot(n,sf_2)
grid on 
hold on 
plot(n, s)



%% functions
function f=fibonacci(n)
% FIBONACCI(n) generates the first n numbers of the fibonacci sequence
f(1)=1;

f(2)=1;
for i=3:n
f(i)=f(i-1)+f(i-2);
end

end