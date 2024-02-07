%checking results
%%far-field formula 5b 
clear all;clc;
warning("off")
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
d1=[100 200];%100:100:400;
d2=100:50:300;

%New params
dl = 100;
ht = 2;%5;%2;
hr = 3;%5;
h = 10;%20;
phi = 0;%pi/2;
shadow_dl = 0.1;


for i=1:length(d1)
    d_far=d1(i);
    Pr = Pt * (Gt * Gr * G * M^2 * N^2 * dx * dy * lambda^2 * F_theta_t *F_theta_r * A_square) ./ (64 * pi^3 * d_far.^2 * d2.^2)* abs((sinc((M * pi / lambda) * (sin(theta_t) * cos(phi_t) + sin(theta_r) * cos(phi_r)) * dx) / sinc((pi / lambda) * (sin(theta_t) * cos(phi_t) + sin(theta_r) * cos(phi_r)) * dx)))^2 * abs((sinc((N * pi / lambda) * (sin(theta_t) * sin(phi_t) + sin(theta_r) * sin(phi_r)) * dy) / sinc((pi / lambda) * (sin(theta_t) * sin(phi_t) + sin(theta_r) * sin(phi_r)) * dy)))^2;
    dbw_power=10*log10(Pr/0.001);

    d = 10 * (2 * h - ht - hr);
    dl = sqrt(d_far.^2 + d2.^2 + 2 .* d_far .* d2 .* cos(theta_t+theta_r));
% %     dl = sqrt((ht - hr)^2 + d^2);
    Pr_dl = sqrt(Gt * Gr) * shadow_dl./ dl;

    Pr_dl_1 = Pt * Gt * Gr * lambda * lambda * shadow_dl * shadow_dl./ (16*pi*pi*(dl).^(2)) ;
    Pr_dl_dB = 10*log10(Pr_dl_1/0.001);%+85

    Pr_ris_fact_2 = (sinc((M * pi / lambda) * (sin(theta_t) * cos(phi_t) + sin(theta_r) * cos(phi_r)) * dx) / sinc((pi / lambda) * (sin(theta_t) * cos(phi_t) + sin(theta_r) * cos(phi_r)) * dx));
    Pr_ris_fact_3 = (sinc((N * pi / lambda) * (sin(theta_t) * sin(phi_t) + sin(theta_r) * sin(phi_r)) * dy) / sinc((pi / lambda) * (sin(theta_t) * sin(phi_t) + sin(theta_r) * sin(phi_r)) * dy));
    Pr_ris_fact_1 = M * N * sqrt(A_square * Gt * Gr * G * F_theta_t * F_theta_r * dx * dy) * exp(1j*phi) ./ (2 * sqrt(pi) * d_far * d2);
% %     Pr_ris_fact_4 = 1-(1j * 4 * pi / (lambda * d)) * (h * h - h * ht - h * hr + ht * hr);
% %     Pr_ris_fact_4 = exp(-1j * (2 * pi / lambda) * (d_far + d2 - dl));
    d_est = 5*(2*h-ht-hr);%sqrt(dl.^2 - (ht-hr).^2);
    Pr_h =  exp(-1j * (2 * pi / lambda) * (2./d_est) * (h * h - h * ht - h * hr + ht * hr));
    Pr_ris_fact_4 = exp(-1j * (2 * pi / lambda) * (2./d_est) * (h * h - h * ht - h * hr + ht * hr));
    Pr_ris = Pt * (lambda^2 / (16 * pi^2)) .* abs(Pr_ris_fact_1 .* Pr_ris_fact_2 .* Pr_ris_fact_3).^2;
    Pr_dl_ris = Pt * (lambda^2 / (16 * pi^2)) .* abs(Pr_dl + (Pr_ris_fact_1 .* Pr_ris_fact_2 .* Pr_ris_fact_3)).^2;
    Pr_ris_h = Pt * (lambda^2 / (16 * pi^2)) .* abs(Pr_ris_fact_1 .* Pr_ris_fact_2 .* Pr_ris_fact_3 .* Pr_ris_fact_4).^2;
    Pr_dl_ris_h = Pt * (lambda^2 / (16 * pi^2)) .* abs(Pr_dl + (Pr_ris_fact_1 .* Pr_ris_fact_2 .* Pr_ris_fact_3 .* Pr_ris_fact_4)).^2;
    dbw_power_1=10*log10(Pr_ris/0.001);
    dbw_power_2=10*log10(Pr_dl_ris/0.001);
    dbw_power_3=10*log10(Pr_ris_h/0.001);
    dbw_power_4=10*log10(Pr_dl_ris_h/0.001);

    figure %(1)
    plot(d2,dbw_power,"Marker","o");
    hold all;
    hold on
    plot(d2,Pr_dl_dB,'-c');
    ylim([-90 -50]);
    ylabel("Received Power in dbm");
    xlabel("When distance between RIS and Rx changes");
    title("Far-Field scenario");
    legend("when d1=100m","when d1=200m","when d1=300m","when d1=400m","Location","best");
    grid on

    figure %(3)
    plot(d2,dbw_power_1,"Marker","o");
    hold all;
    plot(d2,dbw_power_2,"Marker","*");
    plot(d2,dbw_power_3,"Marker","square");
    plot(d2,dbw_power_4,"Marker","pentagram");
    ylim([-90 -50]);
    ylabel("Received Power in dbm");
    xlabel("When distance between RIS and Rx changes");
    title("Far-Field scenario",h);
    legend("when d1=100m","when d1=100m","when d1=200m","when d1=200m","when d1=300m","when d1=300m","when d1=400m","when d1=400m","Location","best");
    grid on

end
%near-field bradcasing case
% % figure (12)
% % d1_near=1:1.5:4;
% % d2_near=20:20:200;
% % for i=1:length(d1_near)
% %     d1=d1_near(i);
% %     Ar=(Gr*lambda_square)/(4*pi);
% %     P_near_received=(Gt * Ar * A_square) ./ (4 * pi * ((d1 + d2_near).^2)) * Pt;
% %     plot(d2_near,10*log10(P_near_received/0.001),"Marker","*");
% %     ylabel("Received Power in dbm");
% %     xlabel("When distance between RIS and Rx changes");
% %     title("Near-Field scenario");
% %     legend("when d1=1m","when d1=2.5m","when d1=4m","Location","best");
% %     hold all;
% %     grid on;
% % end
