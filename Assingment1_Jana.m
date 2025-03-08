%%
% Jana Rusrus
% Assignment 1: estimation of correlation and PSD functions for LTI systems inputs and outputs
%% 
clear all
close all
clc
%
pkg load signal

%% (B.1) Generating the instance of the sequences x( n) , d( n) , y(n)
L=2^16;
g=[0,2,1];
h=[0,0,1,3];
extra_samples_to_remove_transients = max(length(g),length(h))-1;

sigma_x_sq=1;
mx=5;
x=sqrt(sigma_x_sq) * randn(1, L+extra_samples_to_remove_transients)+mx;
d=filter(g,1,x);
y=filter(h,1,x);

% remove transients
x=x(extra_samples_to_remove_transients+1:end);
d=d(extra_samples_to_remove_transients+1:end);
y=y(extra_samples_to_remove_transients+1:end);


%% ------------------ Part 1 (b) ------------------%%

% remove mean
x_hat = x - mean(x);
d_hat = d - mean(d);
y_hat = y - mean(y);


%% ------------------ Phi (yd)------------------%%

figure(1)

%% The theoretical expressions for Phi (yd)
impulse = @(n) n==0;
theoritcal_phi_yd = @(l) 300 + impulse(l) + 5*impulse(l-1) + 6*impulse(l-2);
lag = -10:10; 
subplot(212), stem(lag, theoritcal_phi_yd(lag) );
axis([-10 10 294 308])
grid on
title('Theoritical Cross-correlation \phi_{yd}(l)');
xlabel('lag l');

% Cross -correlation estimate Phi(yd)
[phi_est_yd,lags] = xcorr(y,d,'unbiased');
subplot(211), stem(lags,phi_est_yd)
axis([-10 10 294 308])
grid on
title('Unbiased estimate of Cross-correlation \phi\^_{yd}(l)');
xlabel('lag l');

%% ------------ Phi (yx)------------%%
figure(2)

%% The theoretical expressions for Phi (yx)
impulse = @(n) n==0;
theoritcal_phi_yx = @(l) 100 + impulse(l-2) + 3*impulse(l-3);
lag = -10:10; 
subplot(212), stem(lag, theoritcal_phi_yx(lag) );
axis([-10 10 85 105])
grid on
title('Theoritical Cross-correlation \phi_{yx}(l)');
xlabel('lag l');

% Cross -correlation estimate Phi(yx)
[phi_est_yx,lags] = xcorr(y,x,'unbiased');
subplot(211), stem(lags,phi_est_yx)
axis([-10 10 85 105])
grid on
title('Unbiased estimate of Cross-correlation \phi\^_{yx}(l)');
xlabel('lag l');


%% ------------ Phi (yy)------------%%
figure(3)

%% The theoretical expressions for Phi (yy)
impulse = @(n) n==0;
theoritcal_phi_yy = @(l) 400 + 3*impulse(l+1) + 11*impulse(l) + 3*impulse(l-1);
lag = -10:10; 
subplot(212), stem(lag, theoritcal_phi_yy(lag) );
axis([-10 10 390 415])
grid on
title('Theoritical Auto-correlation \phi_{yy}(l)');
xlabel('lag l');

% Auto -correlation estimate Phi(yy)
[phi_est_yy,lags] = xcorr(y,y,'unbiased');
subplot(211), stem(lags,phi_est_yy)
axis([-10 10 390 415])
grid on
title('Unbiased estimate of Auto-correlation \phi\^_{yy}(l)');
xlabel('lag l');

%% ------------------ Part 1 (c) ------------------%%

%% ------------------ Gamma (yd)------------------%%

figure(4)
%% The theoretical expressions for Gamma (yd)
impulse = @(n) n==0;
theoritcal_gamma_yd = @(l) impulse(l) + 5 * impulse(l-1) + 6*impulse(l-2);
lag = -10:10; 
subplot(212), stem(lag, theoritcal_gamma_yd(lag) );
axis([-10 10 0 10])
grid on
title('Theoritical Cross-covarience \gamma_{yd}(l)');
xlabel('lag l');

% Cross -correlation estimate gamma(yd)
[gamma_est_yd,lags] = xcov(y,d,'biased');
subplot(211), stem(lags,gamma_est_yd)
axis([-10 10 0 10])
grid on
title('biased estimate of Cross-coveriance \gamma\^_{yd}(l)');
xlabel('lag l');

%% ------------------ Gamma (yx)------------------%%
figure(5)
%% The theoretical expressions for Gamma (yx)
impulse = @(n) n==0;
theoritcal_gamma_yx = @(l) impulse(l-2) + 3*impulse(l-3);
lag = -10:10; 
subplot(211), stem(lag, theoritcal_gamma_yx(lag) );
axis([-10 10 0 10])
grid on
title('Theoritical Cross-covarience \gamma_{yx}(l)');
xlabel('lag l');

% Cross -cov estimate gamma(yx)
[gamma_est_yx,lags] = xcov(y,x,'biased');
subplot(212), stem(lags,gamma_est_yx)
axis([-10 10 0 10])
grid on
title('biased estimate of Cross-covariance \gamma\^_{yx}(l)');
xlabel('lag l');


%% ------------------ Gamma (yy)------------------%%
figure(6)
%% The theoretical expressions for Gamma (yy)
impulse = @(n) n==0;
theoritcal_gamma_yy = @(l)  3 * impulse(l+1) + 10*impulse(l) + 3 * impulse(l-1);
lag = -10:10; 
subplot(212), stem(lag, theoritcal_gamma_yy(lag) );
axis([-10 10 -2 12])
grid on
title('Theoritical Auto-covarience \gamma_{yy}(l)');
xlabel('lag l');

% Cross -cov estimate gamma(yy)
[gamma_est_yy,lags] = xcov(y, 'biased');
subplot(211), stem(lags,gamma_est_yy)
axis([-10 10 -2 12])
grid on
title('biased estimate of auto-covariance \gamma\^_{yy}(l)');
xlabel('lag l');


%% ------------------ Part 1 (d) ------------------%%
L = 2^16;
sigma_x_sq = 1;
mx = 5;
% Generate 100 sequence y(n)
x100 = sqrt(sigma_x_sq) * randn(100, L+extra_samples_to_remove_transients) + mx;
y100 = filter(h,1,x100')';

x100 = x100(:, extra_samples_to_remove_transients+1:end);
y100 = y100(:, extra_samples_to_remove_transients+1:end);

phi_est_yy100 = [];
gamma_est_yy100 = [];

for i = 1:100
  [phi, lags_phi] = xcorr(y100(i,:) ,'unbiased');
  phi_est_yy100 = [phi_est_yy100; phi];

  [gamma, lags_gamma] = xcov(y100(i,:), 'biased');  
  gamma_est_yy100 = [gamma_est_yy100; gamma];

end

mu_phi_est = mean(phi_est_yy100);
var_phi_est = mean( (phi_est_yy100 - repmat(mu_phi_est, 100, 1)).^2 );

mu_gamma_est = mean(gamma_est_yy100);
var_gamma_est = mean( (gamma_est_yy100 - repmat(mu_gamma_est, 100, 1)).^2 );

figure(7)
subplot(211), plot(lags_phi, var_phi_est)
grid on
title('Varience of Unbiased estimate of Auto-correlation \phi\^_{yy}(l)');
xlabel('lag l');
%axis([round(-2*pi*1e4) round(2*pi*1e4) 0 5])
ylim([0,5])

subplot(212), plot(lags_gamma, var_gamma_est)
grid on
title('biased estimate of auto-covariance \gamma\^_{yy}(l)');
xlabel('lag l');
%axis([round(-2*pi*1e4) round(2*pi*1e4) 0 4e-3])
ylim([0, 4e-3])

%% ------------------ Part 2 ------------------%%
L = [12, 16];

figure(8)

for i = 1:2
    % Generate 100 sequence y(n)
    x100 = sqrt(sigma_x_sq) * randn(100, 2^L(i)+extra_samples_to_remove_transients) + mx;
    d100 = filter(g,1,x100')';
    y100 = filter(h,1,x100')';

    x100 = x100(:, extra_samples_to_remove_transients+1:end);
    d100 = d100(:, extra_samples_to_remove_transients+1:end);
    y100 = y100(:, extra_samples_to_remove_transients+1:end);

    
    % Periodogram estimate
    Pyy_100 = [];
    for k = 1:100
      [gamma, lag] = xcov(y100(k,:), 'biased');
      
      fft_len = max(4*length(gamma), 4096);
      w = linspace(-pi, pi, fft_len);
      
      Gamma_w = fft( hanning(length(gamma))' .* gamma, ...
                     fft_len ...
                   );                
      Gamma_w = fftshift( abs(Gamma_w) );
      
      Pyy_100 = [Pyy_100; Gamma_w];
      
    end

    % Variance
    Pyy_mean = mean(Pyy_100);
    Pyy_var = mean( (Pyy_100 - repmat(Pyy_mean,100,1)).^2 );

    subplot(2,1,i), plot(w, Pyy_var)
    xlabel('\omega [rad/sample]')
    title(sprintf('Variance of Periodogram Auto-PSD Estimate P_{yy} , L = 2^{%d}', L(i)))
    grid on

end

%% ------------------ Part 3 ------------------%%

L = 2^16;
M = 128;
k = 1023;
overlap = 0.5* M;

fft_len = max(4*M, 4096); 
w = linspace(-pi, pi, fft_len);

norm_window = hamming(M)' / sqrt(mean(hamming(M).^2));

%% yd

% Ideal signal
s1 =  (exp(2*j*w) + 3*exp(j*w)) .* (2*exp(-j*w) + exp(-2*j*w));   % Phi_yd

% Estimate 
Pyd_est = [];

for i = 1 : M - overlap : length(y) - overlap
  y_seg = y_hat( i : i+M-1 ) .* norm_window;
  Y_seg = fft( y_seg, fft_len);     % FFT for the noisy signal
  Y_seg = abs(fftshift(Y_seg));
  
  d_seg = d_hat( i : i+M-1 ) .* norm_window;
  D_seg = fft( d_seg , fft_len);     % FFT for the noisy signal
  D_seg = abs(fftshift(D_seg));
  
  Gamma_seg = Y_seg .* conj(D_seg) / M;
  Pyd_est = [ Pyd_est ; Gamma_seg ];
  
end

Pyd_welch = mean(Pyd_est);

figure(9)   
subplot(2,1,1), plot(w,abs(s1))
xlabel('\omega [rad/sample]')
title('|cross PSD \Phi_{yd}|')
grid on

subplot(2,1,2), plot(w, abs(Pyd_welch))
xlabel('\omega [rad/sample]')
title('Welch Estimate P_{yd}')
grid on


%% yx
% Ideal signal
s2 = exp(j*2* w)+ 3 * exp(j*3*w);

% Estimate
Pyx_est = [];

for i = 1 : M - overlap : length(y) - overlap
  y_seg = y_hat( i : i+M-1 ) .* norm_window;
  Y_seg = fft( y_seg, fft_len);     % FFT for the noisy signal
  Y_seg = abs(fftshift(Y_seg));
  
  x_seg = x_hat( i : i+M-1 ) .* norm_window;
  X_seg = fft( x_seg , fft_len);     % FFT for the noisy signal
  X_seg = abs(fftshift(X_seg));
  
  Gamma_seg = Y_seg .* conj(X_seg) / M;
  Pyx_est = [ Pyx_est ; Gamma_seg ];
  
end

Pyx_welch = mean(Pyx_est);

figure(10)    
subplot(2,1,1), plot(w,abs(s2))
xlabel('\omega [rad/sample]')
title('|cross PSD \Phi_{yx}|')
grid on

subplot(2,1,2), plot(w, abs(Pyx_welch))
xlabel('\omega [rad/sample]')
title('Welch Estimate P{yx}')
grid on


%% yy
% Ideal signal
s3 = 6* cos(w)+ 10;
%Estimate
Pyy_est = [];

for i = 1 : M - overlap : length(y) - overlap
  y_seg = y_hat( i : i+M-1 ) .* norm_window;
  Y_seg = fft( y_seg, fft_len);     % FFT for the noisy signal
  Y_seg = abs(fftshift(Y_seg));
    
  Gamma_seg = Y_seg .* conj(Y_seg) / M;
  Pyy_est = [ Pyy_est ; Gamma_seg ];
  
end

Pyy_welch = mean(Pyy_est);

figure(11)    
subplot(2,1,1), plot(w,abs(s3))
xlabel('\omega [rad/sample]')
title('|Auto PSD \Phi_{yy}|')
grid on

subplot(2,1,2), plot(w, abs(Pyy_welch))
xlabel('\omega [rad/sample]')
title('Welch Estimate P{yy}')
grid on


%% Pxx
Pxx_est = [];

for i = 1 : M - overlap : length(y) - overlap
  x_seg = x_hat( i : i+M-1 ) .* norm_window;
  X_seg = fft( x_seg, fft_len);     % FFT for the noisy signal
  X_seg = abs(fftshift(X_seg));
    
  Gamma_seg = X_seg .* conj(X_seg) / M;
  Pxx_est = [ Pxx_est ; Gamma_seg ];
  
end

Pxx_welch = mean(Pxx_est);

%% Coherence

Coherence = Pyx_welch ./ sqrt(Pyy_welch .* Pxx_welch);

figure(12)      
plot(w, abs(Coherence))
xlabel('\omega [rad/sample]')
title('Coherence \Gamma_{yx}')
grid on


