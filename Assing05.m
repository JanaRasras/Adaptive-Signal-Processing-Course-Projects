%%
% Jana Rusrus
% Assignment 5: NLMS, AP, LS and RLS algorithms

%% 
clear all
close all
clc
%
% pkg load signal

tic

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


u =@(l) l>=0;

Weiner.phi_xx =@(l) 25.2525   * (0.98.^l .* u(l) + 0.98.^(-l).* u(-l-1));
Weiner.phi_dx =@(l) 824.9166 * 0.98.^l - 799.6641414 * 0.95.^l .* u(l) + ...
                    365.978626 * 0.98.^(-l).* u(-l-1) + 340.7261 * 0.95 .^l .* u(l);

P = Weiner.phi_dx(0 : N-1)';
R = toeplitz(Weiner.phi_xx(0: N-1)');


%% W opt
wopt = R\P;

% Eigenvalue Spread
eigen_max = max(eig(R));
eigen_min = min(eig(R));
spread = eigen_max/eigen_min;

%%  Parameters 

nb_iter = 10000;

% Smoothing (2nd Plot)
LPF = ones(1,N)/N;     % short: measure speed of convergence

% NLMS
NLMS.mu_Range = linspace(0, 2, 10);
%NLMS.mu = 1.6;                               % Q.1 Largest conv
%NLMS.mu = 1.0;                               % Q.1 Largest conv
NLMS.mu = 1.0;                               % Q.3,4,5  adjus

NLMS.MSE = zeros(nb_iter, 1);                 % Figure 1
NLMS.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
NLMS.learning_curve = zeros(nb_iter, 1);      % Figure 3

NLMS.w = zeros(N,1);

% AP
M = 5; 
AP.X = zeros(N, M);                         % Augmented-X
AP.d = zeros(1, M);                         % Augmented-d

AP.mu_Range = linspace(0, 2, 10);
AP.mu = 0.29;                               % Q.3,Q5 adjust
%AP.mu = 1.9;                               % Q1. largest conv
%AP.mu = NLMS.mu;                             %Q.4

AP.MSE = zeros(nb_iter, 1);
AP.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
AP.learning_curve = zeros(nb_iter, 1);      % Figure 3

AP.w = zeros(N,1);


% LS
N = 100;
LS.delta = 1e+2;                          % adjust/ Forgetting factor
%LS.lambda = 0.999;                        % Q3,4
LS.lambda = 0.9855;                          %Q5

LS.MSE = zeros(nb_iter, 1);                 % Figure 1
LS.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
LS.learning_curve = zeros(nb_iter, 1);      % Figure 3

LS.psi = eye(N) * LS.delta;
LS.phi = zeros(N,1);
LS.w = zeros(N,1);


% RLS
RLS.lambda = LS.lambda;                     %Q.4
RLS.lambda = 0.9855;                           %Q.5
RLS.delta = 1/ LS.delta;
RLS.w = zeros(N, 1);
RLS.k = zeros(N,1);
RLS.inv_R = eye(N) * RLS.delta;

RLS.MSE = zeros(nb_iter, 1);                 % Figure 1
RLS.SmoothMSE = zeros(nb_iter, 1);           % Figure 2
RLS.learning_curve = zeros(nb_iter, 1);      % Figure 3

%% Iterative

k = 10 * N;
for n = k : k+nb_iter-1 
    
    %---- NLMS -----------------%
    err1 = dA(n) - NLMS.w' * xA(n:-1:n-N+1);
    err3 = d_power - NLMS.w'*P - NLMS.w'*conj(P) + NLMS.w' * R * NLMS.w;
    
    NLMS.MSE(n-k+1) = err1.^2 / d_power;                     
    NLMS.learning_curve(n-k+1) = abs(err3) / d_power;
    
    %update
    NLMS.w = NLMS.w + NLMS.mu / ( 1e-10 + xA(n:-1:n-N+1)' * xA(n:-1:n-N+1) ) * xA(n:-1:n-N+1) * err1'; 
    
    %---- AP -----------------%
    AP.X(:,2:M) = AP.X(:,1:M-1);                % Shift-left
    AP.X(:,1) = [xA(n); AP.X(1:end-1,1)];       % Prepend new x (end=N)
    
    AP.d(2:M) = AP.d(1:M-1);
    AP.d(1) = dA(n);
    
    err1 = AP.d - AP.w' * AP.X;
    err3 = d_power - AP.w'*P - AP.w'*conj(P) + AP.w' * R * AP.w;
    
    AP.MSE(n-k+1) = err1(1).^2 / d_power;                     
    AP.learning_curve(n-k+1) = abs(err3) / d_power;
    
    %update
    dw = AP.X * inv(AP.X' * AP.X + eye(M)* 1e-10 )* err1';
    AP.w = AP.w + AP.mu * dw; 
    
    %if n == 5000
    %    AP.mu = AP.mu / 100;
    %end
    
    %---- LS -----------------% 
    err1 = dA(n) - LS.w' * xA(n:-1:n-N+1);
    err3 = d_power - LS.w'*P - LS.w'*conj(P) + LS.w' * R * LS.w;
    
    LS.psi = LS.lambda * LS.psi + xA(n:-1:n-N+1) *xA(n:-1:n-N+1)';
    LS.phi = LS.lambda * LS.phi + xA(n:-1:n-N+1)* dA(n)';
    LS.psi_inv = inv(LS.psi);
    LS.w = LS.psi_inv * LS.phi;
    
    LS.MSE(n-k+1) = err1.^2 / d_power;                     
    LS.learning_curve(n-k+1) = abs(err3) / d_power;
    
    %---- RLS -----------------% 
    err1 = dA(n) - RLS.w' * xA(n:-1:n-N+1);
    err3 = d_power - RLS.w'*P - RLS.w'*conj(P) + RLS.w' * R * RLS.w;
    
    RLS.u = RLS.inv_R * xA(n:-1:n-N+1);
    RLS.k = RLS.u / (RLS.lambda  + xA(n:-1:n-N+1)' * RLS.u);
    RLS.inv_R = 1/ RLS.lambda * (RLS.inv_R - RLS.k * RLS.u');
    up = triu(RLS.inv_R);
    RLS.inv_R = up + up';
    RLS.inv_R =RLS.inv_R - diag( 0.5 * real( diag(RLS.inv_R)) + sqrt(-1) * imag( diag(RLS.inv_R) )); % pg30 l140
    
    RLS.w = RLS.w + RLS.k * err1';
    
    RLS.MSE(n-k+1) = err1.^2 / d_power;                     
    RLS.learning_curve(n-k+1) = abs(err3) / d_power;
    
    % Reference (Weiner)
    err4 = d_power - (conj(P)' * (R \ P));
    err4_learning_curve(n-k+1) = err4 / d_power;
end

err2 = filter(LPF, 1, NLMS.MSE); 
NLMS.SmoothMSE = err2; 


err2 = filter(LPF, 1, AP.MSE); 
AP.SmoothMSE = err2;

err2 = filter(LPF, 1, LS.MSE); 
LS.SmoothMSE = err2;


err2 = filter(LPF, 1, RLS.MSE); 
RLS.SmoothMSE = err2;

misadj1 = abs(10*log10(err4_learning_curve(end)) - 10*log10(NLMS.SmoothMSE(end)));
misadj2 = abs(10*log10(err4_learning_curve(end)) - 10*log10(AP.SmoothMSE(end)));
misadj3 = abs(10*log10(err4_learning_curve(end)) - 10*log10(LS.SmoothMSE(end)));
misadj4 = abs(10*log10(err4_learning_curve(end)) - 10*log10(RLS.SmoothMSE(end)));

%%
figure(1)
subplot(311); plot(1:nb_iter, 10*log10(NLMS.MSE),...
                   1:nb_iter,10*log10(AP.MSE), ...
                   1:nb_iter,10*log10(LS.MSE+eps),...
                   1:nb_iter,10*log10(RLS.MSE+eps), '--',...
                   1:nb_iter, 10*log10(err4_learning_curve), 'k');
               
 
               
S_NLMS = sprintf('NLMS, mu norm=%0.2f ', NLMS.mu );
S_AP = sprintf('AP, mu norm=%0.2f ', AP.mu );
S_LS = sprintf('LS');
S_RLS = sprintf('RLS');
S_MMSE = sprintf('MMSE');

legend(S_NLMS, S_AP, S_LS, S_RLS,S_MMSE);

ylim([-105.0 5.0])
%xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized squared error (non-smoothed)')
grid on


indx = 1:nb_iter-N+6;
subplot(312); plot( indx, 10*log10(NLMS.SmoothMSE(N-5:end)'),...
                    indx, 10*log10(AP.SmoothMSE(N-5:end))', ...
                    indx, 10*log10(LS.SmoothMSE(N-5:end)+eps)', ...
                    indx, 10*log10(RLS.SmoothMSE(N-5:end)+eps)', '--',...
                    1:nb_iter, 10*log10(err4_learning_curve), 'k');
                    

S_NLMS = sprintf('NLMS, mu norm=%0.2f ', NLMS.mu );
S_AP = sprintf('AP, mu norm=%0.2f ', AP.mu );
S_LS = sprintf('LS');
S_RLS = sprintf('RLS');

legend(S_NLMS, S_AP, S_LS,S_RLS, S_MMSE);


ylim([-45.0 5.0])
%xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized squared error (smoothed)')
grid on

subplot(313); plot( 1:nb_iter, 10*log10(NLMS.learning_curve),...
                   1:nb_iter, 10*log10(AP.learning_curve), ...
                   1:nb_iter, 10*log10(LS.learning_curve+eps),...
                   1:nb_iter, 10*log10(RLS.learning_curve+eps), '--',...
                   1:nb_iter, 10*log10(err4_learning_curve), 'k');

S_NLMS = sprintf('NLMS, mu norm=%0.2f', NLMS.mu );
S_AP = sprintf('AP, mu norm=%0.2f ', AP.mu );
S_LS = sprintf('LS');
S_RLS = sprintf('RLS');

legend(S_NLMS, S_AP, S_LS, S_RLS, S_MMSE);


ylim([-45.0 5.0])
%xlim([0 2000])
xlabel('iterations')
ylabel('dB')
title('Normalized MSE learning curve')
grid on

toc


