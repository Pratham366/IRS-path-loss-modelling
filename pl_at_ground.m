%checking results
%%far-field formula 5b 
clear all;clc;
Gt=125;
Gr=125;
G=1;
Pt=0.01;
M=102;
N=100;
dx=0.01;
dy=0.01;
theta_t=pi/4;
theta_r=pi/4;
phi_t=pi;
phi_r=0;
f=10.5e12;
lambda=0.0286;
lambda_square=(0.0286)^2;
F_theta_t=cos(theta_t)^3;
F_theta_r=cos(theta_r)^3;
A_square=0.9^2;
d1=100:100:400;
d2=100:50:300;
for i=1:length(d1)
    d_far=d1(i);
    Pr = Pt * (Gt * Gr * G * M^2 * N^2 * dx * dy * lambda^2 * F_theta_t *F_theta_r * A_square) ./ (64 * pi^3 * d_far.^2 * d2.^2)* abs((sinc((M * pi / lambda) * (sin(theta_t) * cos(phi_t) + sin(theta_r) * cos(phi_r)) * dx) / sinc((pi / lambda) * (sin(theta_t) * cos(phi_t) + sin(theta_r) * cos(phi_r)) * dx)))* abs((sinc((N * pi / lambda) * (sin(theta_t) * sin(phi_t) + sin(theta_r) * sin(phi_r)) * dy) / sinc((pi / lambda) * (sin(theta_t) * sin(phi_t) + sin(theta_r) * sin(phi_r)) * dy)))^2;
    dbw_power=10*log10(Pr/0.001);
    figure (1)
    plot(d2,dbw_power,"Marker","o");
    hold all;
    ylim([-90 -50]);
    ylabel("Received Power in dbm");
    xlabel("When distance between RIS and Rx changes");
    title("Far-Field scenario");
    legend("when d1=100m","when d1=200m","when d1=300m","when d1=400m","Location","best");
    grid on
end
%near-field bradcasing case
figure (2)
d1_near=1:1.5:4;
d2_near=20:20:200;
for i=1:length(d1_near)
    d1=d1_near(i);
    Ar=(Gr*lambda_square)/(4*pi);
    P_near_received=(Gt * Ar * A_square) ./ (4 * pi * ((d1 + d2_near).^2)) * Pt;
    plot(d2_near,10*log10(P_near_received/0.001),"Marker","*");
    ylabel("Received Power in dbm");
    xlabel("When distance between RIS and Rx changes");
    title("Near-Field scenario");
    legend("when d1=1m","when d1=2.5m","when d1=4m","Location","best");
    hold all;
    grid on;
end
