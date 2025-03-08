%%
% Jana Rusrus
% Assignment 3: Source separation,
%% 
clear
close all
clc

% pkg load signal

%% Question 1

% reading signals from sources A, b and microphones 1, 2
[data.mic1.a, fs] = audioread('mic1_SA_only_interval1.wav');
[data.mic1.b, fs] = audioread('mic1_SB_only_interval2.wav');
[data.mic1.ab, fs] = audioread('mic1_SA_and_SB_interval3.wav');

[data.mic2.a, fs] = audioread('mic2_SA_only_interval1.wav');
[data.mic2.b, fs] = audioread('mic2_SB_only_interval2.wav');
[data.mic2.ab, fs] = audioread('mic2_SA_and_SB_interval3.wav');


%% ------------------ Source A ------------------%
N = 10;
lag = length(data.mic1.a);

% desired signal
SA.dn = data.mic1.a;  
% P Matrix
phi_x1d = xcorr(data.mic1.a, SA.dn, 'biased');
phi_x2d = xcorr(data.mic2.a, SA.dn, 'biased');

% R Matrix
phi_x1x1 = xcorr(data.mic1.ab, data.mic1.ab, 'biased');                
phi_x1x2 = xcorr(data.mic1.ab, data.mic2.ab, 'biased');                
phi_x2x1 = xcorr(data.mic2.ab, data.mic1.ab, 'biased');                
phi_x2x2 = xcorr(data.mic2.ab, data.mic2.ab, 'biased'); 

% Wiener
Q1.SA.P = [phi_x1d(lag:-1:lag-N+1); phi_x2d(lag:-1:lag-N+1)];

subR11 = toeplitz( phi_x1x1(lag: -1: lag-N+1) , phi_x1x1(lag: lag+N-1) );
subR12 = toeplitz( phi_x1x2(lag: -1: lag-N+1) , phi_x1x2(lag: lag+N-1) );
subR21 = toeplitz( phi_x2x1(lag: -1: lag-N+1) , phi_x2x1(lag: lag+N-1) );
subR22 = toeplitz( phi_x2x2(lag: -1: lag-N+1) , phi_x2x2(lag: lag+N-1) );
R = [subR11 subR12 ; subR21 subR22];

Q1.SA.Wopt = R \ Q1.SA.P;

figure(); stem(Q1.SA.Wopt, 'filled', 'MarkerSize',3)
title('Wopt-SA'); 

% filtered output
Q1.SA.output = filter(Q1.SA.Wopt(1:N), 1, data.mic1.ab) + filter(Q1.SA.Wopt(N+1:end), 1, data.mic2.ab);
%sound(Q1.SA.output, fs)
audiowrite('Q1_SA_Output.wav', 0.9 * Q1.SA.output/max(abs(Q1.SA.output)), fs);

%% ------------------ Source B ------------------%
% desired signal
SB.dn = data.mic1.b;  

phi_x1d = xcorr(data.mic1.b, SB.dn, 'biased');
phi_x2d = xcorr(data.mic2.b, SB.dn, 'biased');

% Wiener
Q1.SB.P = [phi_x1d(lag:-1:lag-N+1); phi_x2d(lag:-1:lag-N+1)];

Q1.SB.Wopt = R \ Q1.SB.P;

figure(); stem(Q1.SB.Wopt, 'filled', 'MarkerSize',3)
title('Wopt-SB'); 

% filtered output
Q1.SB.output = filter(Q1.SB.Wopt(1:N), 1, data.mic1.ab) + filter(Q1.SB.Wopt(N+1:end), 1, data.mic2.ab);
sound(Q1.SB.output, fs)
audiowrite('Q1_SB_Output.wav',Q1.SB.output,fs);


%% Power Ratio
% Source A
xa = filter(Q1.SA.Wopt(1:N), 1, data.mic1.a) + filter(Q1.SA.Wopt(N+1:end), 1, data.mic2.a);                    
xb = filter(Q1.SA.Wopt(1:N), 1, data.mic1.b) + filter(Q1.SA.Wopt(N+1:end), 1, data.mic2.b);                    
%SNR.A = 10*log10((xa'*xa)/(xb'*xb));

ga1 = (xa' * xa) / ((data.mic1.a)' * (data.mic1.a));
ga2 = (xb' * xb) / ((data.mic1.b)' * (data.mic1.b));

Q1.SIR.A = 10*log10(ga1/ga2);

% Source B
xa = filter(Q1.SB.Wopt(1:N), 1, data.mic1.a) + filter(Q1.SB.Wopt(N+1:end), 1, data.mic2.a);                    
xb = filter(Q1.SB.Wopt(1:N), 1, data.mic1.b) + filter(Q1.SB.Wopt(N+1:end), 1, data.mic2.b);                    
%SNR.B = 10*log10((xb'*xb)/(xa'*xa));

gb1 = (xa' * xa) / ((data.mic1.a)' * (data.mic1.a));
gb2 = (xb' * xb) / ((data.mic1.b)' * (data.mic1.b));

Q1.SIR.B = 10*log10(gb2/gb1);

Q1.SIR
%% Question 2
SNR = 30; % note: if SNR is very large (ex. >100 dB), it is as if there is no noise (Q1)
%SNR = 20;
%SNR = 10;
%SNR = 0;

rng(13); % set random generator seed, so that the noise sequences remain the same between
         % different runs, allowing a more direct comparison

speechA_mic1_interval1_rms = sqrt(sum(abs(data.mic1.a).^2) / length(data.mic1.a));
noise = randn(length(data.mic1.a), 1);
noise = noise*speechA_mic1_interval1_rms / (10.^(SNR/20));
data.mic1.a_noisy = data.mic1.a + noise; 

speechA_mic2_interval1_rms = sqrt(sum(abs(data.mic2.a).^2) / length(data.mic2.a));
noise = randn(length(data.mic2.a),1);
noise = noise*speechA_mic2_interval1_rms / (10.^(SNR/20));
data.mic2.a_noisy = data.mic2.a + noise;

speechB_mic1_interval2_rms = sqrt(sum(abs(data.mic1.b).^2)/length(data.mic1.b));
noise = randn(length(data.mic1.b),1);
noise = noise*speechB_mic1_interval2_rms/(10.^(SNR/20));
data.mic1.b_noisy = data.mic1.b + noise;

speechB_mic2_interval2_rms = sqrt(sum(abs(data.mic2.b).^2) / length(data.mic2.b));
noise = randn(length(data.mic2.b),1);
noise = noise*speechB_mic2_interval2_rms / (10.^(SNR/20));
data.mic2.b_noisy = data.mic2.b + noise;

speechAB_mic1_interval3_rms = sqrt(sum(abs(data.mic1.ab).^2)/length(data.mic1.ab));
noise=randn(length(data.mic1.ab),1);
noise=noise*speechAB_mic1_interval3_rms / (10.^(SNR/20));
data.mic1.ab_noisy = data.mic1.ab + noise;

speechAB_mic2_interval3_rms = sqrt(sum(abs(data.mic2.ab).^2) / length(data.mic2.ab));
noise=randn(length(data.mic2.ab),1);
noise=noise*speechAB_mic2_interval3_rms / (10.^(SNR/20));
data.mic2.ab_noisy = data.mic2.ab+noise;

%% ------------------ Source A ------------------%

% desired signal
SA.dn = data.mic1.a;  

phi_x1d = xcorr(data.mic1.a_noisy, SA.dn, 'biased');
phi_x2d = xcorr(data.mic2.a_noisy, SA.dn, 'biased');

phi_x1x1 = xcorr(data.mic1.ab_noisy, data.mic1.ab_noisy, 'biased');                
phi_x1x2 = xcorr(data.mic1.ab_noisy, data.mic2.ab_noisy, 'biased');                
phi_x2x1 = xcorr(data.mic2.ab_noisy, data.mic1.ab_noisy, 'biased');                
phi_x2x2 = xcorr(data.mic2.ab_noisy, data.mic2.ab_noisy, 'biased'); 

% Wiener
Q2.SA.P = [phi_x1d(lag:-1:lag-N+1); phi_x2d(lag:-1:lag-N+1)];

subR11 = toeplitz( phi_x1x1(lag: -1: lag-N+1) , phi_x1x1(lag: lag+N-1) );
subR12 = toeplitz( phi_x1x2(lag: -1: lag-N+1) , phi_x1x2(lag: lag+N-1) );
subR21 = toeplitz( phi_x2x1(lag: -1: lag-N+1) , phi_x2x1(lag: lag+N-1) );
subR22 = toeplitz( phi_x2x2(lag: -1: lag-N+1) , phi_x2x2(lag: lag+N-1) );
R = [subR11 subR12 ; subR21 subR22];

Q2.SA.Wopt = R \ Q2.SA.P;

% filtered output
Q2.SA.output = filter(Q2.SA.Wopt(1:N), 1, data.mic1.ab_noisy) + filter(Q2.SA.Wopt(N+1:end), 1, data.mic2.ab_noisy);
%sound(SA.output, fs)
audiowrite('Q2_SA_Output.wav', 0.9 * Q2.SA.output/max(abs(Q2.SA.output)), fs);

%% ------------------ Source B ------------------%
% desired signal
SB.dn = data.mic1.b;  

phi_x1d = xcorr(data.mic1.b_noisy, SB.dn, 'biased');
phi_x2d = xcorr(data.mic2.b_noisy, SB.dn, 'biased');

% Wiener
Q2.SB.P = [phi_x1d(lag:-1:lag-N+1); phi_x2d(lag:-1:lag-N+1)];

Q2.SB.Wopt = R \ Q2.SB.P;

% filtered output
Q2.SB.output = filter(Q2.SB.Wopt(1:N), 1, data.mic1.ab_noisy) + filter(Q2.SB.Wopt(N+1:end), 1, data.mic2.ab_noisy);
%sound(Q2.SB.output, fs)
audiowrite('Q2_SB_Output.wav',Q2.SB.output, fs);


%% Power Ratio
% Source A
xa = filter(Q2.SA.Wopt(1:N), 1, data.mic1.a_noisy) + filter(Q2.SA.Wopt(N+1:end), 1, data.mic2.a_noisy);                    
xb = filter(Q2.SA.Wopt(1:N), 1, data.mic1.b_noisy) + filter(Q2.SA.Wopt(N+1:end), 1, data.mic2.b_noisy);                    

ga1 = (xa' * xa) / ((data.mic1.a_noisy )' * (data.mic1.a_noisy ));
ga2 = (xb' * xb) / ((data.mic1.b_noisy )' * (data.mic1.b_noisy ));

Q2.SIR.A = 10*log10(ga1/ga2);

% Source B
xa = filter(Q2.SB.Wopt(1:N), 1, data.mic1.a_noisy) + filter(Q2.SB.Wopt(N+1:end), 1, data.mic2.a_noisy);                    
xb = filter(Q2.SB.Wopt(1:N), 1, data.mic1.b_noisy) + filter(Q2.SB.Wopt(N+1:end), 1, data.mic2.b_noisy);                    

gb1 = (xa' * xa) / ((data.mic1.a_noisy )' * (data.mic1.a_noisy ));
gb2 = (xb' * xb) / ((data.mic1.b_noisy )' * (data.mic1.b_noisy ));

Q2.SIR.B = 10 * log10(gb2/gb1);

Q2.SIR