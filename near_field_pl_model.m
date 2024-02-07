
%%far-field


%%near-field
N = 100;
M = 102;
dx = 0.01;
dy = 0.01;
A = 0.9;
lambda = 0.0286;
theta_t=pi/4;
theta_r=pi/4;

% % d_1 = 2;
d_2 = [1:0.5:5];
d1_near = 2; 
d2_near=20:10:200;

Gt=125;
Gr=125;
Pt=0.01;
f=10.5e12;
lambda=0.0286;
lambda_square=(0.0286)^2;

Pr = Pt * Gt * Gr * A * A * lambda * lambda ./ (16*pi*pi*((d1_near+d2_near).^(2)));
Pr_dB = 10*log10(Pr/0.001);%+85

%%direct link Near field broadcasting case
ht = 2;
hr = 1;
h_arr = [20 15 10 5];
shadow_dl = 1;
dl = sqrt(d1_near^(2)+d2_near.^(2)++ 2 .* d1_near .* d2_near .* cos(theta_t+theta_r));
d_est = 5*(2*h-ht-hr);%sqrt((d_1+d_2).^2 - (2*h-ht-hr).^2)+4;
% % dl = sqrt((ht-hr)*(ht-hr)+d_est.^2);
phi = pi/2;

Pr_dl = Pt * Gt * Gr * lambda * lambda ./ (16*pi*pi*(dl).^(2)) ;
Pr_dl_dB = 10*log10(Pr_dl/0.001);%+85

Pr_comb = Pt * Gt * Gr * ((lambda * lambda ./(16*pi*pi))) .* (abs((shadow_dl./dl) + ((A*exp(1j*0)./(d1_near+d2_near))))).^(2);
Pr_comb_dB = 10*log10(Pr_comb/0.001);%+85;

Pr_comb_1 = Pt * Gt * Gr * ((lambda * lambda ./(16*pi*pi))) .* (abs((shadow_dl./dl) + ((A*exp(1j*pi/4)./(d1_near+d2_near))))).^(2);
Pr_comb_1_dB = 10*log10(Pr_comb_1/0.001);%+85;

Pr_comb_2 = Pt * Gt * Gr * ((lambda * lambda ./(16*pi*pi))) .* (abs((shadow_dl./dl) + ((A*exp(1j*pi/2)./(d1_near+d2_near))))).^(2);
Pr_comb_2_dB = 10*log10(Pr_comb_2/0.001);%+85;

Pr_comb_3 = Pt * Gt * Gr * ((lambda * lambda ./(16*pi*pi))) .* (abs((shadow_dl./dl) + ((A*exp(1j*3*pi/4)./(d1_near+d2_near))))).^(2);
Pr_comb_3_dB = 10*log10(Pr_comb_3/0.001);%+85;

% % Pr_comb_2 = Pt * Gt * Gr * [(lambda * lambda ./(16*pi*pi))] .* (abs([1./dl] + [A./(d1+d2)].*(1-1j*4*(pi/(lambda*d))*[h*h+ht*hr-h*ht-h*hr]))).^(2);
% % Pr_comb_2_dB = 10*log(Pr_comb_2);%+85;

figure
grid on
plot(d2_near,Pr_dB,'-r');
hold on
plot(d2_near,Pr_dl_dB,'-c');
hold on
plot(d2_near,Pr_comb_dB,'-b')
hold on
plot(d2_near,Pr_comb_1_dB,'-g')
hold on
plot(d2_near,Pr_comb_2_dB,'-k')
hold on
plot(d2_near,Pr_comb_3_dB,'-m')
xlabel('d2')
ylabel('Received Power (in dBm)')
legend("ris","dl","\phi=0","\phi=pi/4","\phi=pi/2","\phi=3pi/4","Location","best");
title("Near-Field scenario");
grid on

% % Pr_h =  1-1j*((4*pi ./(lambda*d_est)))*(h*h+ht*hr-h*ht-h*hr);
% % dl = 1;
for i=1:length(h_arr)
h=h_arr(i);
d_est = 5*(2*h-ht-hr);%sqrt(dl.^2 - (ht-hr).^2);
Pr_h =  exp(-1j * (2 * pi / lambda) * (2./d_est) * (h * h - h * ht - h * hr + ht * hr));
Pr_comb = Pt * Gt * Gr * [(lambda * lambda ./(16*pi*pi))] .* (abs([shadow_dl./dl] + [(A*exp(1j*0)./(d1_near+d2_near)).*Pr_h])).^(2);
Pr_comb_h_dB = 10*log10(Pr_comb/0.001);%+85;
Pr_comb_1 = Pt * Gt * Gr * [(lambda * lambda ./(16*pi*pi))] .* (abs([shadow_dl./dl] + [(A*exp(1j*pi/4)./(d1_near+d2_near)).*Pr_h])).^(2);
Pr_comb_h_1_dB = 10*log10(Pr_comb_1/0.001);%+85;
Pr_comb_2 = Pt * Gt * Gr * [(lambda * lambda ./(16*pi*pi))] .* (abs([shadow_dl./dl] + [(A*exp(1j*pi/2)./(d1_near+d2_near)).*Pr_h])).^(2);
Pr_comb_h_2_dB = 10*log10(Pr_comb_2/0.001);%+85;
Pr_comb_3 = Pt * Gt * Gr * [(lambda * lambda ./(16*pi*pi))] .* (abs([shadow_dl./dl] + [(A*exp(1j*3*pi/4)./(d1_near+d2_near)).*Pr_h])).^(2);
Pr_comb_h_3_dB = 10*log10(Pr_comb_3/0.001);%+85;

figure
grid on
plot(d2_near,Pr_dB,'-r');
hold on
plot(d2_near,Pr_dl_dB,'-c');
hold on
plot(d2_near,Pr_comb_h_dB,'-b')
hold on
plot(d2_near,Pr_comb_h_1_dB,'-g')
hold on
plot(d2_near,Pr_comb_h_2_dB,'-k')
hold on
plot(d2_near,Pr_comb_h_3_dB,'-m')
xlabel('d2')
ylabel('Received Power (in dBm)')
legend("ris","dl","\phi=0","\phi=pi/4","\phi=pi/2","\phi=3pi/4","Location","best");
title("Near-Field scenario with height",h);
grid on
end
