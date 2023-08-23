
%% DFT Leakage Minimization (filtering) using Hamming window function
clc;    clear;   close all;
   Fs = 64000;                  % Sampling frequency (24KHz)
   ts = 1/Fs;     % seconds
   DFT_points = 64;      N = DFT_points;
     a = 8;   Fc_1 = 3300;   b = 6;   Fc_2 = 3700;
     ind = 1;   x = [];
    for n = 1:DFT_points
        m = n-1;
        x1(ind) = a*sin(2*pi*Fc_1*m*ts);
        x2(ind) = b*sin(2*pi*Fc_2*m*ts);
        ind = ind + 1;
    end
    x_comb = x1 + x2;     
    t = 1:DFT_points;
% Plot first N discrete values of x_comb signal:
   figure(1);   plot(t,x_comb,'b--o');   grid on;
   xlabel('Time (millisecond)');     ylabel('Signal amplitude')
   title('x_comb signal versus time');   zoom xon;
   
%% DFT of x_comb (using exponential equation):
Dft_x_comb = zeros(N,1);     
Dft_x_comb = dft(x_comb, N);

mf = 0:DFT_points-1;
figure(2); 
stem(mf,Dft_x_comb,'LineStyle','--',...
     'MarkerSize',15,'Marker','s',...
     'MarkerFaceColor','black',...
     'MarkerEdgeColor','green')
grid on;
title('Magnitude of Dft_x_comb')
xlabel('m (KHz)')
ylabel('Magnitude')

% Filtering using Hamming Window
w_Ham = zeros(N,1);
for ii = 1:DFT_points
    n = ii - 1;
    w_Ham(ii,1) = 0.54 - 0.46*cos(2*pi*(n/(N-1)));      %Define Hamming window function
end

figure(3);
stem(t,w_Ham); grid on;
title('Hamming Window');  xlabel('n');  ylabel('w')
x_sig = x_comb';
x_Ham = x_sig.*w_Ham;               % Multiplication of Hamming Window and sin function: x_comb

figure(4)
stem(t,x_Ham); grid on;
title('Multiplication of Hamming Window and sin function');
xlabel('Time');  ylabel('Amplitude')
x_sig_Ham = x_Ham';
 X_dft_Ham = dft(x_sig_Ham, N);     % DFT of x_sig_Ham

 figure(5); 
stem(mf,X_dft_Ham,'LineStyle','--',...
     'MarkerSize',15,'Marker','s',...
     'MarkerFaceColor','black',...
     'MarkerEdgeColor','green')
grid on;
title('Magnitude of X_dft_Ham')
xlabel('m (KHz)')
ylabel('Magnitude')

   