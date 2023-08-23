
%% A DFT shifting theorem example (using exponential form of DFT equation)
clc;   clear;   close all;

DFT_points = 8;      j = sqrt(-1);
N = DFT_points;   % The number of samples of input sequence- x and the number of frequency points in the DFT output  
k = 3;            % Amount of shift in Time-domain signal  
Fs = 8e3;         % Sampling frequency 8KHz
ts = 1/Fs;
ind = 1;   x = [];
for n = 1:N+k
    m = n-1;
    x(ind) = sin(2*pi*1000*m*ts)+0.5*sin(2*pi*2000*m*ts+(3*pi/4));
    ind = ind + 1;
end

t = 1:N+k;
% Plot first N discrete values of x_seq1 signal:
   figure(1);   plot(t,x,'b--o');   grid on;
   xlabel('Time (millisecond)');     ylabel('Signal amplitude')
   title('x signal versus time');   zoom xon;

x_1 = x(1,1:N); 

x_2 = x(1,1+k:end);            % time sequence shifted by k samples
   
f_analysis_1m = (1*Fs)/N;      % DFT analysis frequency for m=1

X_exp = zeros(N,1);
X_exp = dft(x_1, N);

% for m = 1:DFT_points
%     X_1 = zeros(N,1);   m_e = m - 1;
%     for ii = 1:N
%         n = ii-1;   X_p1 = [];   
%         X_p1 = x_1(1,ii)*exp((-j*2*pi*n*m_e)/N);                    
%         X_1(ii,1) = X_p1;  
%     end
%     X_exp(m,1) = sum(X_1(:,1));
% end

X_shifted = zeros(N,1);     
X_shifted = dft(x_2, N);

% for m = 1:DFT_points
%     X_1 = zeros(N,1);   m_e = m - 1;
%     for ii = 1:N
%         n = ii-1;   X_p1 = [];   
%         X_p1 = x_2(1,ii)*exp((-j*2*pi*n*m_e)/N);                    
%         X_1(ii,1) = X_p1;  
%     end
%     X_shifted(m,1) = sum(X_1(:,1));
% end

% m = 1;
X_dft_k = zeros(N,1);
Xe_real = zeros(N,1);     Xe_imag = zeros(N,1);
for m = 1:N
    m_e = m-1;
    X_dft_k(m,1) = exp((j*2*pi*k*m_e)/N)*X_exp(m,1);
    
    Xe_real(m,1) = real(X_dft_k(m,1));
            if Xe_real(m,1) > 0 && Xe_real(m,1) < 1e-10
               Xe_real(m,1) = 0;  
            end  
            if Xe_real(m,1) < 0 && Xe_real(m,1) > -1e-10
               Xe_real(m,1) = 0;  
            end
        Xe_imag(m,1) = imag(X_dft_k(m,1));
            if Xe_imag(m,1) > 0 && Xe_imag(m,1) < 1e-10
               Xe_imag(m,1) = 0;  
            end  
            if Xe_imag(m,1) < 0 && Xe_imag(m,1) > -1e-10
               Xe_imag(m,1) = 0;  
            end
        X_dft_k_mag(m,1) = sqrt(Xe_real(m,1).^2 + Xe_imag(m,1).^2);
       
end


mf = 0:DFT_points-1;
figure(4); 
stem(mf,X_shifted,'LineStyle','--',...
     'MarkerSize',15,'Marker','s',...
     'MarkerFaceColor','black',...
     'MarkerEdgeColor','green')
grid on;
title('Magnitude of X_shifted(m)')
xlabel('m (KHz)')
ylabel('Magnitude')

figure(5); 
stem(mf,X_dft_k_mag,'LineStyle','--',...
     'MarkerSize',15,'Marker','s',...
     'MarkerFaceColor','blue',...
     'MarkerEdgeColor','green')
grid on;
title('Magnitude of X_dft_k(m)')
xlabel('m (KHz)')
ylabel('Magnitude')



% Verify DFT shifting property: 
DFT_Linearity_error = max(abs(X_shifted - X_dft_k_mag));
    if DFT_Linearity_error < 1e-9
       'DFT_shifting_property_is_proved'  
    end



















