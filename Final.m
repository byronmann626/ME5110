close all;clear;clc;
%Givens
k1 = 30000; k2 = 30000;
c1 = 3000; c2 = 3000;
L1 = 1; L2 = 1.5;
lam = 5;
v0 = 50/3600*1000;
x0 = [0 0];
A = 0.01;
m = 2000;
Ic = 2500;
T = lam/v0;

%Road Profile
r1 = @(t) A*sin(2*pi/T*t);
r2 = @(t) A*sin(2*pi/T*t+pi);
r1_dot = @(t) A*(2*pi/T)*cos(2*pi/T*t);
r2_dot = @(t) A*(2*pi/T)*cos(2*pi/T*t+pi);

%K,C, A, M Matrices
M =[m 0; 0 Ic];
K = [k1+k2 k1*L1-k2*L2; k1*L1-k2*L2 k1*L1^2+k2*L2^2];
C = [c1+c2 c1*L1-c2*L2; c1*L1-c2*L2 c1*L1^2+c2*L2^2];
A = [zeros(2,2) eye(2); -M\K -M\C];

%Question 2 & 3
% Undamped Nat Frequency (Cholesky)
L = chol(M,'Lower');
ktilde = inv(L)*K*inv(L');
[Un,wn] = eig(ktilde);
Us = inv(L')*Un;
Un = [Us(:,1)/Us(2,1) Us(:,2)/Us(2,2)];

%Eigen Vectors and Mode Shapes (Damped and Undamped)
[V,D] = eig(A);
omega = diag(D);
U = V(1:2,:);
U = real([U(:,1)/U(2,1) U(:,3)/U(2,3)]);
tgrid = 0:0.01:5; tgrid = tgrid(:);

%Question 4:
B = []; C1 = eye(4); D = [];
sys = ss(A,B,C1,D);
%Mode Shape 1 Response
s0_sys1 = [U(:,1);0;0];
sys1_resp = initial(sys,s0_sys1,tgrid);
%Mode Shape 2 Response
s0_sys2 = [U(:,2);0;0];
sys2_resp = initial(sys,s0_sys2,tgrid);
%Plot
figure('Name','Problem 4: Individual Mode Shape Response','NumberTitle','off');
subplot(2,1,1)
yyaxis left
plot(tgrid,sys1_resp(:,1));
ylabel('Bounce Displacement (m)');
yyaxis right
plot(tgrid,sys1_resp(:,2));
ylabel('Pitch (rad)');
xlabel('Time (s)');
title('First Mode Shape');
subplot(2,1,2)
yyaxis left
plot(tgrid,sys2_resp(:,1));
ylabel('Displacement (m)');
yyaxis right
plot(tgrid,sys2_resp(:,2));
ylabel('Pitch (rad)');
xlabel('Time (s)');
title('Second Mode Shape');

%Question 5: 
%Free Response for y(0) = 0.014 theta(0) = 0.05
s0_free = [0.014 0.05 0 0]';
free_response = initial(sys, s0_free, tgrid);
figure('Name','Problem 5: Free Response for y(0) = 0.014 theta(0) = 0.05','NumberTitle','off');
subplot(2,1,1)
yyaxis left
plot(tgrid,free_response(:,1));
ylabel('Bounce Displacement (m)');
yyaxis right
plot(tgrid,free_response(:,2));
ylabel('Pitch Displacement (rad)');
legend('Mode Shape 1', 'Mode Shape 2');
xlabel('Time(s)');
subplot(2,1,2)
yyaxis left
plot(tgrid,free_response(:,3))
ylabel('Bounce velocity (m/s)');
yyaxis right
plot(tgrid,free_response(:,4));
ylabel('Pitch velocity (rad/s)');
legend('Mode Shape 1', 'Mode Shape 2');
xlabel('Time(s)');

% %Question 6: 
s = tf('s');

T = [m*s^2 + (c1+c2)*s + (k1 + k2) , (c1*L1 - c2*L2)*s + (k1*L1 - k2*L2);
    (L1*c1 - L2*c2)*s + (L1*k1 - L2*k2) , Ic*s^2 + ((L2^2) *c2 + (L1^2) *c1)*s + ((L2^2)*k2 + (L1^2) *k1)];

H = [(k1 + c1*s) (k2 + c2*s);
    (k1*L1 + c1*L1*s) (k2*L2 - c2*L2*s)];


R = ((m*s^2 + (c1+c2)*s + (k1 + k2))*(Ic*s^2 + ((L2^2) *c2 + (L^2) *c1)*s + ((L2^2) *k2 + (L1^2) *k1)) - ((L1*c1 - L2*c2)*s + (L1*k1 - L2*k2))*((c1*L1 - c2*L2)*s + (k1*L1 - k2*L2)));

Tmod = [Ic*s^2 + ((L2^2) *c2 + (L1^2) *c1)*s + ((L2^2) *k2 + (L1^2) *k1) -((c1*L1 - c2*L2)*s + (k1*L1 - k2*L2));
    -((L1*c1 - L2*c2)*s + (L1*k1 - L2*k2))  m*s^2 + (c1+c2)*s + (k1 + k2)];

TT1 = T\H



fprintf('\n\nProblem 6:\n\n');
fprintf('Transfer Functions:\n')
fprintf('Y(s)/R1(s):\n')
TT1(1,1)
fprintf('Y(s)/R2(s):\n')
TT1(1,2)
fprintf('Theta(s)/R1(s):\n')
TT1(2,1)
fprintf('Theta(s)/R2(s):\n')
TT1(2,2)

nf = 500;
fmin = .01;
fmax = 1000;
fvec = logspace(log10(fmin), log10(fmax),nf);
wvec = 2*pi*fvec;

for m = 1:2
    
    sys = TT1(1,m);
    [mag,phase] = bode(sys,wvec);
    mag = squeeze(mag); phase = squeeze(phase);
    % max(mag)
    % max(phase)
    
    figure('name','Bounce Output Responses')
    subplot(211), semilogx(fvec,20*log10(mag), 'Linewidth',2)
    xlabel('Frequency [Hz]')
    ylabel('Mag (db of m/m)')
    title(sprintf('Frequency Response Magnitude: Y(s)/R%i(s)',m))
    grid
    
    subplot(212), semilogx(fvec,phase,'Linewidth',2)
    xlabel('Frequency [Hz]')
    ylabel('Phase [deg]')
    title(sprintf('Frequency Response Phase: Y(s)/R%i(s)',m))
    grid
end


for n = 1:2
    sys2 = TT1(2,n);
    [mag,phase] = bode(sys2,wvec);
    mag = squeeze(mag); phase = squeeze(phase);
    
    figure('name','Pitch Output Responses')
    subplot(211), semilogx(fvec,20*log10(mag), 'Linewidth',2)
    xlabel('Frequency [Hz]')
    ylabel('Mag (db of m/m)')
    title(sprintf('Frequency Response Magnitude: Theta(s)/R%i(s)',n))
    grid
    
    subplot(212), semilogx(fvec,phase,'Linewidth',2)
    xlabel('Frequency [Hz]')
    ylabel('Phase [deg]')
    title(sprintf('Frequency Response Phase: Theta(s)/R%i(s)',n))
    grid
end


%Question 7: 
%Numerical Simulation of Response (forced)
tforce = 0:0.01:10; tforce = tforce(:);
k_force = [k1 k2; -k1*L1 k2*L2];
c_force = [c1 c2; -c1*L1 c2*L2];
s0 = zeros(4,1);
F = @(t) [k_force*[r1(t);r2(t)]+c_force*[r1_dot(t);r2_dot(t)]];
f = @(t,s) [s(3); s(4); M\(F(t)-K*[s(1);s(2)]-C*[s(3);s(4)])];
[t,s] = ode45(f,tforce,s0);
figure('Name','Problem 7: Forced Vibration','NumberTitle','off');
subplot(2,1,1)
yyaxis left
plot(t,s(:,1));
ylabel('Bounce Displacement (m)');
yyaxis right
plot(t,s(:,2));
ylabel('Pitch Displacement (m)');
legend('Mode Shape 1', 'Mode Shape 2');
xlabel('Time(s)');
subplot(2,1,2)
yyaxis left
plot(t,s(:,3));
ylabel('Bounce Velocity (m/s)');
yyaxis right
plot(t,s(:,4));
ylabel('Pitch Velocity (rad/s)');
legend('Mode Shape 1', 'Mode Shape 2');
xlabel('Time(s)');


%Question 8:


