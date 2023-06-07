%Valsa Code
%% Author: Kurt Cantrell 
%% Last updated: 5/23/22
% The purpose of this script is to generate R and C values for any order
% with a Valsa FOE circuit topology.

clear all;
close all;
clc;

%%% Parameters %%%

%
alpha = 0.5; %input('What is the desired alpha for the expression? ');
m = 10; %input('What is the desired order for the expression? '); % Order Default 1/3

Rin = 10000; %input('What is the desired R1 for the expression? '); %R value for time constant Default 10k
Cin = 0.000001; %input('What is the desired C1 for the expression? '); %C value for time constant Default 1 u

deltPhi = 0.5; %input('What is the allowable ripple for the expression? ');; %Allowable ripple Default 0.5
Dp = 1e4; %Default 1e4


t = Rin*Cin;


f = logspace(-1,9, 4001);
w = f*2*pi;
s = 1i*w;


j = 1i;

phiDegrees = -90*alpha;
phiRadians = -alpha*pi/2;




ab = 0.24/(1+ deltPhi);
a = 10^(alpha*log10(ab));
%a = inv(a);
b = ab/a;

Rp = (Rin*(1-a))/(a);
Cp = (Cin*b^m)/(1-b);

for n=1:1:m
    R(n) = Rin*a^(n-1);

    C(n) = Cin*b^(n-1);
end

k = ceil(m/2);

w_av = ((a/b)^0.25)*(1/(Rin*Cin*(ab)^2));

x=0;
for n = 1:1:m
    x = x + (j*w_av*C(n))/(1+j*w_av*R(n)*C(n));
end

Yjw_av = 1/Rp + j*w_av*Cp + x;

Z_av = 1/abs(Yjw_av);

D = Z_av*(w_av^(alpha));

temp = Dp/D;

Rp = Rp*temp;
Cp = Cp/temp;

for n=1:1:m
    R(n) = R(n)*temp;

    C(n) = C(n)/temp;
end


Ztot = (Rp*(1./(s*Cp)))./(Rp+(1./(s*Cp)));


for n = 1:m

    Ztot = (Ztot.*(R(n) + 1./(s*C(n))))./(Ztot+(R(n) + 1./(s*C(n))));

end



%{
figure (3) %error
subplot(2,1,1)
semilogx((w/(2*pi)), Zmag);
xlabel('Frequency (Hz)'); ylabel('Error (%)');
title('FI: Magnitude Error for Alpha of', num2str(a));
ylim([-5 5]);
%

subplot(2,1,2)
semilogx((w/(2*pi)), Zang);
xlabel('Frequency (Hz)'); ylabel('Error (%)');
title('FI: Phase Error for Alpha of', num2str(a));
ylim([-5 5]);%% plotting verification

%}
figure (1) %magnitdue 
%loglog((w/(2*pi)), (abs(Zs))); hold on;
%loglog((w/(2*pi)), (abs(Zab))); hold on;
loglog((w/(2*pi)), (abs(Ztot))); %%this fixed it. log 10 
%legend('Ideal', 'A & B Coefficients', 'R and C Values');
xlabel('Frequency (Hz)'); ylabel('Magnitude (Ohm)');
title(['FI: Impedance Mag with Alpha of ', num2str(alpha), ' and Capacitance of ', num2str(a),'F']);
grid on;

%
figure (2) %phase 
%semilogx((w/(2*pi)), (angle(Zs)*180/pi)); hold on;
%semilogx((w/(2*pi)), (angle(Zab)*180/pi)); hold on;
semilogx((w/(2*pi)), (angle(Ztot)*180/pi));
%legend('Ideal', 'A & B Coefficients', 'R and C Values');
xlabel('Frequency (Hz)'); ylabel('Phase Angle (degrees)');
title(['FI: Impedance Phase with Alpha of ', num2str(alpha), ' and Capacitance of ', num2str(a),'F']);
grid on;
%}

R = R';
C = C';