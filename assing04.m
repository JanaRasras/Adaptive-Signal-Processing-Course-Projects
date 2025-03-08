%%
% Jana Rusrus
% Assignment 4: Steepest descent, Newton, LMS and Newton-LMS algorithm

%% 
clear all
close all
clc
%
% pkg load signal

%% Generating the signals
rng(13); % to use same random sequence between different simulations
L = 10000000;
u = randn(L+200,1); % white Gaussian, variance of one
x = filter(1,[1 -0.98],u); % t(n) filter
x = x(101:end); % discard first 100 samples to remove transients
% in simulated signals
v = randn(L+100,1); % white Gaussian, variance of one
dA = filter(1,[1 -0.95],x) + v; % h(n) filter plus additive noise
dA = dA(101:end); % discard first 100 samples to remove transients
% in simulated signals
xA = x(101:end); % discard first 100 samples to remove transients
% in simulated signals


%% Finding R and p
N = 100;            % coefficients of the filter
d_power =  7249.25; %dA(N: N + nb_iter)' * dA(N: N + nb_iter)/ nb_iter;

%{
[phi_est_xx,lagsx] = xcorr(xA, xA,'biased');
[phi_est_dx,lagsd] = xcorr(dA, xA,'biased');
[phi_est_dd,lagsdd] = xcorr(dA, dA,'biased');

% W opt 
n1 = find(lagsx == 0);
n2 = find(lagsd == 0);
n3 = find(lagsdd == 0);

P = phi_est_dx(n2:n2+N-1);
R = toeplitz(phi_est_xx(n1:n1+N-1));
%}


u =@(l) l>=0;

Weiner.phi_xx =@(l) 25.2525   * (0.98.^l .* u(l) + 0.98.^(-l).* u(-l-1));
Weiner.phi_dx =@(l) 824.9166 * 0.98.^l - 799.6641414 * 0.95.^l .* u(l) + ...
                    365.978626 * 0.98.^(-l).* u(-l-1) + 340.7261 * 0.95 .^l .* u(l);

P = Weiner.phi_dx(0 : N-1)';
R = toeplitz(Weiner.phi_xx(0: N-1)');


%{
[phi_est_xx,lagsx] = xcorr(xA, xA,'unbiased');
[phi_est_dx,lagsd] = xcorr(dA, xA,'unbiased');
[phi_est_dd,lagsdd] = xcorr(dA, dA,'unbiased');

n1 = find(lagsx == 0);
n2 = find(lagsd == 0);
n3 = find(lagsdd == 0);

P2 = phi_est_dx(n2:n2+N-1);
R2 = toeplitz(phi_est_xx(n1:n1+N-1));
%}

%% W opt

%{
load p
load R
P = p;
%}

wopt = R\P;

% Eigenvalue Spread
eigen_max = max(eig(R));
eigen_min = min(eig(R));
spread = eigen_max/eigen_min;

%% Parameters
nb_iter = 3000;

% SD
SD.mu_Range = linspace(0, 2/eigen_max, 10);
SD.mu = 0.001; %SD.mu_Range(9);                     % adjust

SD.MSE = zeros(nb_iter, 1);                 % Figure 1
SD.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
SD.learning_curve = zeros(nb_iter, 1);      % Figure 3

SD.w = zeros(N,1);

% Newton
Newton.mu_max = 2;
Newton.mu = 1;            % adjust

Newton.MSE = zeros(nb_iter, 1);                 % Figure 1
Newton.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
Newton.learning_curve = zeros(nb_iter, 1);      % Figure 3

Newton.w = zeros(N,1);


% LMS
LMS.mu_Range = linspace(0, 2 /(3* trace(R)), 10);
LMS.mu = 0.5/10000 *2; %LMS.mu_Range(3); 

LMS.MSE = zeros(nb_iter, 1);                 % Figure 1
LMS.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
LMS.learning_curve = zeros(nb_iter, 1);      % Figure 3

LMS.w = zeros(N,1);

% Newton_LMS
Newton_LMS.mu_Range = linspace(0, 2/N, 10);
Newton_LMS.mu = 12/1000 ;            % adjust

Newton_LMS.MSE = zeros(nb_iter, 1);                 % Figure 1
Newton_LMS.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
Newton_LMS.learning_curve = zeros(nb_iter, 1);      % Figure 3

Newton_LMS.w = zeros(N,1);

% Smoothing (2nd Plot)
LPF = ones(1,N)/N;     % short: measure speed of convergence

%% Iterative
k = 10*N;

for n = k : k+nb_iter-1          % here k = n i.e. one iteration per sample
    
    %---- Newton -----------------
    % Evaluate
    err1 = dA(n) - Newton.w' * xA(n:-1:n-N+1);
    err3 = d_power - Newton.w'*P - Newton.w'*conj(P) + Newton.w' * R * Newton.w;
    
    Newton.MSE(n-k+1) = err1 ^ 2 / d_power;
    Newton.learning_curve(n-k+1) = abs(err3) / d_power;
   
    % Update
    dw = R \ (R * Newton.w - P);
    Newton.w = Newton.w - Newton.mu * dw;
    
    %----  SD  -----------------   
    % Evaluate
    err1 = dA(n) - SD.w' * xA(n:-1:n-N+1);
    err3 = d_power - conj(SD.w)'*P - SD.w'*conj(P) + conj(SD.w)' * R * SD.w;
    
    SD.MSE(n-k+1) = err1 ^ 2 / d_power;
    SD.learning_curve(n-k+1) = abs(err3) / d_power;
    
    % Update
    dw = SD.mu * ( R * SD.w - P );
    SD.w = SD.w -  dw;
    
    %---- LMS -----------------
    % Evaluate
    err1 = dA(n) - LMS.w' * xA(n:-1:n-N+1);
    err3 = d_power - LMS.w'*P - LMS.w'*conj(P) + LMS.w' * R * LMS.w;
    
    LMS.MSE(n-k+1) = err1 ^ 2 / d_power;
    LMS.learning_curve(n-k+1) = abs(err3) / d_power;
    
    % Update
    dw = xA(n:-1:n-N+1) * conj(err1);
    LMS.w = LMS.w + LMS.mu * dw;
    
    %---- Newton_LMS -----------------
    % Evaluate
    err1 = dA(n) - Newton_LMS.w' * xA(n:-1:n-N+1);
    err3 = abs( d_power - Newton_LMS.w'*P - Newton_LMS.w'*conj(P) + Newton_LMS.w' * R * Newton_LMS.w);
    
    Newton_LMS.MSE(n-k+1) = err1 ^ 2 / d_power;
    Newton_LMS.learning_curve(n-k+1) = abs(err3) / d_power;

    % Update
    dw = R \ (xA(n:-1:n-N+1) * conj(err1)');
    Newton_LMS.w = Newton_LMS.w + Newton_LMS.mu * dw;
    
   
    
    %---- Plot W_est -----------------
    %{
    figure(1)
    subplot(321); stem(wopt), ylim([0 1.2])
    title('W opt');
    
    subplot(323); stem(SD.w), ylim([0 1.2])
    title('SD');
    
    subplot(324); stem(Newton.w), ylim([0 1.2])
    title('Newton');
    
    subplot(325); stem(LMS.w), ylim([0 1.2])
    title('LMS');
    
    subplot(326); stem(Newton_LMS.w), ylim([0 1.2])
    title('Newton LMS');
    pause(0.001)
    %}
end
% SD
err2 = filter(LPF, 1, SD.MSE); 
SD.SmoothMSE = err2; %(length(LPF):end);     % ignore 1st L-1 samples

% Newton
err2 = filter(LPF, 1, Newton.MSE);
Newton.SmoothMSE = err2; % (length(LPF):end);     % ignore 1st L-1 samples

% LMS
err2 = filter(LPF, 1, LMS.MSE); 
LMS.SmoothMSE = err2; % (length(LPF):end);     % ignore 1st L-1 samples

% Newton_LMS
err2 = filter(LPF, 1, Newton_LMS.MSE); 
Newton_LMS.SmoothMSE = err2; %(length(LPF):end);     % ignore 1st L-1 samples



%% Plot
figure(2)

subplot(311); plot(1:nb_iter, 10*log10(SD.MSE), ...
                   1:nb_iter, 10*log10(LMS.MSE),...
                   1:nb_iter, 10*log10(Newton.MSE), ...
                   1:nb_iter, 10*log10(Newton_LMS.MSE))
                
S_SD = sprintf('SD, mu norm=%0.2f e-3', SD.mu * 1000);
s_Newton = sprintf('Newton, mu norm=%0.2f', Newton.mu);
S_LMS = sprintf('LMS, mu norm=%0.2f e-3', LMS.mu * 1000);
s_NLMS = sprintf('Newton-LMS, mu norm=%0.2f e-3', Newton_LMS.mu * 1000);

legend(S_SD, S_LMS, s_Newton, s_NLMS);

ylim([-105.0 5.0])
xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized squared error (non-smoothed)')
grid on

indx = 1:nb_iter-N+6;
subplot(312); plot( indx, 10*log10(SD.SmoothMSE(N-5:end)), ...
                    indx, 10*log10(LMS.SmoothMSE(N-5:end)),...
                    indx, 10*log10(Newton.SmoothMSE(N-5:end) ), ...
                    indx, 10*log10(Newton_LMS.SmoothMSE(N-5:end)))

legend(S_SD,S_LMS, s_Newton,  s_NLMS);

ylim([-45.0 5.0])
xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized squared error (smoothed)')
grid on

subplot(313); plot( 1:nb_iter, 10*log10(SD.learning_curve),...
                    1:nb_iter, 10*log10(LMS.learning_curve),...
                    1:nb_iter, 10*log10(Newton.learning_curve ), ...
                    1:nb_iter, 10*log10(Newton_LMS.learning_curve))

legend(S_SD, S_LMS, s_Newton,  s_NLMS);

ylim([-45.0 5.0])
xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized MSE learning curve')
grid on


%% Questin 3
nb_iter = 5000;

LMS.mu = 0.5/10000 *2;                         % Best case
LMS.MSE = zeros(nb_iter, 1);                 % Figure 1
LMS.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
LMS.learning_curve = zeros(nb_iter, 1);      % Figure 3
LMS.w = zeros(N,1);

SD.mu = LMS.mu;                             % Copy LMS
SD.MSE = zeros(nb_iter, 1);                 % Figure 1
SD.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
SD.learning_curve = zeros(nb_iter, 1);      % Figure 3
SD.w = zeros(N,1);

for n = k : k+nb_iter-1 
    %----  SD  -----------------   
    % Evaluate
    err1 = dA(n) - SD.w' * xA(n:-1:n-N+1);
    err3 = d_power - conj(SD.w)'*P - SD.w'*conj(P) + conj(SD.w)' * R * SD.w;
    
    SD.MSE(n-k+1) = err1 ^ 2 / d_power;
    SD.learning_curve(n-k+1) = abs(err3) / d_power;
    
    % Update
    dw = SD.mu * ( R * SD.w - P );
    SD.w = SD.w -  dw;
    
    %---- LMS -----------------
    % Evaluate
    err1 = dA(n) - LMS.w' * xA(n:-1:n-N+1);
    err3 = d_power - LMS.w'*P - LMS.w'*conj(P) + LMS.w' * R * LMS.w;
    
    LMS.MSE(n-k+1) = err1 ^ 2 / d_power;
    LMS.learning_curve(n-k+1) = abs(err3) / d_power;
    
    % Update
    dw = xA(n:-1:n-N+1) * conj(err1);
    LMS.w = LMS.w + LMS.mu * dw;
    
    % Reference (Weiner)
    err4 = d_power - (conj(P)' * (R \ P));
    err4_learning_curve(n-k+1) = err4 / d_power;
end

% Smooth Error
% SD
err2 = filter(LPF, 1, SD.MSE); 
SD.SmoothMSE = err2; %(length(LPF):end);     % ignore 1st L-1 samples

% LMS
err2 = filter(LPF, 1, LMS.MSE); 
LMS.SmoothMSE = err2; % (length(LPF):end);     % ignore 1st L-1 samples

%% Plot
figure(3)

subplot(311); plot(1:nb_iter, 10*log10(SD.MSE), ...
                   1:nb_iter, 10*log10(LMS.MSE),...
                   1:nb_iter, 10*log10(err4_learning_curve), 'k')
                
S_SD = sprintf('SD');
S_LMS = sprintf('LMS');
S_MMSE = sprintf('MMSE');

legend(S_SD, S_LMS, S_MMSE);

ylim([-105.0 5.0])
%xlim([0 5000])
xlabel('iterations')
ylabel('dB')
title('Normalized squared error (non-smoothed)')
grid on

indx = 1:nb_iter-N+6;
subplot(312); plot( indx, 10*log10(SD.SmoothMSE(N-5:end)), ...
                    indx, 10*log10(LMS.SmoothMSE(N-5:end)), ...
                    1:nb_iter, 10*log10(err4_learning_curve), 'k')

legend(S_SD, S_LMS, S_MMSE);

ylim([-45.0 5.0])
%xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized squared error (smoothed)')
grid on

subplot(313); plot( 1:nb_iter, 10*log10(SD.learning_curve),...
                    1:nb_iter, 10*log10(LMS.learning_curve),...
                    1:nb_iter, 10*log10(err4_learning_curve), 'k')

legend(S_SD, S_LMS, S_MMSE);

ylim([-45.0 5.0])
%xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized MSE learning curve')
grid on

%
misadj = abs(10*log10(err4_learning_curve(end)) - 10*log10(LMS.learning_curve(end)));


%% Questin 4
nb_iter = 10000;

LMS.mu = 0.5/10000 *2;                         % Best case
LMS.MSE = zeros(nb_iter, 1);                 % Figure 1
LMS.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
LMS.learning_curve = zeros(nb_iter, 1);      % Figure 3
LMS.w = zeros(N,1);

SD.mu = LMS.mu;                             % Copy LMS
SD.MSE = zeros(nb_iter, 1);                 % Figure 1
SD.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
SD.learning_curve = zeros(nb_iter, 1);      % Figure 3
SD.w = zeros(N,1);

for n = k : k+nb_iter-1 
    if n-k == 2200
        LMS.mu = LMS.mu / 10;
        SD.mu = SD.mu / 10;
    end    

    %----  SD  -----------------   
    % Evaluate
    err1 = dA(n) - SD.w' * xA(n:-1:n-N+1);
    err3 = d_power - conj(SD.w)'*P - SD.w'*conj(P) + conj(SD.w)' * R * SD.w;
    err4 = d_power - (conj(P)' * (R \ P));
    
    SD.MSE(n-k+1) = err1 ^ 2 / d_power;
    SD.learning_curve(n-k+1) = abs(err3) / d_power;
    
    % Update
    dw = SD.mu * ( R * SD.w - P );
    SD.w = SD.w -  dw;
    
    %---- LMS -----------------
    % Evaluate
    err1 = dA(n) - LMS.w' * xA(n:-1:n-N+1);
    err3 = d_power - LMS.w'*P - LMS.w'*conj(P) + LMS.w' * R * LMS.w;
    
    LMS.MSE(n-k+1) = err1 ^ 2 / d_power;
    LMS.learning_curve(n-k+1) = abs(err3) / d_power;
    
    % Update
    dw = xA(n:-1:n-N+1) * conj(err1);
    LMS.w = LMS.w + LMS.mu * dw;
        
    %---- MMSE -----------------
    err4_learning_curve(n-k+1) = err4 / d_power;
    
end

% Smooth Error
% SD
err2 = filter(LPF, 1, SD.MSE); 
SD.SmoothMSE = err2; %(length(LPF):end);     % ignore 1st L-1 samples

% LMS
err2 = filter(LPF, 1, LMS.MSE); 
LMS.SmoothMSE = err2; % (length(LPF):end);     % ignore 1st L-1 samples

%% Plot
figure(4)

subplot(311); plot(1:nb_iter, 10*log10(SD.MSE), ...
                   1:nb_iter, 10*log10(LMS.MSE),...
                   1:nb_iter, 10*log10(err4_learning_curve), 'k')
                
S_SD = sprintf('SD');
S_LMS = sprintf('LMS');
S_MMSE = sprintf('MMSE');

legend(S_SD, S_LMS, S_MMSE);

ylim([-105.0 5.0])
%xlim([0 5000])
xlabel('iterations')
ylabel('dB')
title('Normalized squared error (non-smoothed)')
grid on

indx = 1:nb_iter-N+6;
subplot(312); plot( indx, 10*log10(SD.SmoothMSE(N-5:end)), ...
                    indx, 10*log10(LMS.SmoothMSE(N-5:end)), ...
                    1:nb_iter, 10*log10(err4_learning_curve), 'k')

legend(S_SD, S_LMS, S_MMSE);

ylim([-45.0 5.0])
%xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized squared error (smoothed)')
grid on

subplot(313); plot( 1:nb_iter, 10*log10(SD.learning_curve),...
                    1:nb_iter, 10*log10(LMS.learning_curve),...
                    1:nb_iter, 10*log10(err4_learning_curve), 'k')

legend(S_SD, S_LMS, S_MMSE);

ylim([-45.0 5.0])
%xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized MSE learning curve')
grid on

% 
misadj = 10*log10(err4_learning_curve(end)) - 10*log10(LMS.learning_curve(end));
