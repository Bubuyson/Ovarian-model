clear;clc;close all
h = 0.01;
t = 0:h:90;
%%
a = 0.6; b = 0.004; c1 = 0.0045; c2 = 0.077; c3 = 0.006; c4 = 0.008;
c5 = 0.045; d1 = 0.5; d2 = 0.8; k1 = 0.6; k2 = 0.5; k3 = 0.6; k4 = 0.65;
%%
e0 = 48; e1 = 0.7; e2 = 2.1; e3 = 1.7; p1 = 0.55; p2 = 0.45; h0 = 270;
h1 = 2.5; h2 = 2; h3 = 10; h4 = 14;
%%
FSH(1) = 142.165;
LH(1) = 28.007;
%%
   MsF(1) = 0;
   GrF(1) = 0;
   PrF(1) = 0;
   Sc1(1) = 0;
   Sc2(1) = 0; 
   Lut1(1) = 0;
   Lut2(1) = 0;
   Lut3(1) = 0;   
   Lut4(1) = 0;
   E2(1) = 50;
   P4(1) = 1;
   Ih(1) = 290;
%    MsF(i) = MsF(i-1) + h * dMsF;
%    GrF(i) = GrF(i-1) + h * dGrF;
%    PrF(i) = PrF(i-1) + h*dPrF(i);
%    Sc1(i) = Scl(i-1) + h*dScl;
%    Sc2(i) = Sc2(i-1) + h*dSc2; 
%    Lut1(i) = Lut1(i-1) + h*dLut1;
%    Lut2(i) = Lut2(i-1) + h*dLut2;
%    Lut3(i) = Lut3(i-1) + h*dLut3;   
%    Lut4(i) = Lut4(i-1) + h*dLut4;

%%
for i = 2:length(t)
   val = 175*exp(-(t(i)-5)^2/120) + 150*exp(-(t(i)-35)^2/160);
   FSH(i) = 250 - (250*(t(i) - 15)^4 )/(1+(t(i)-15)^4) + val;
   LH(i) = 380 - 352*(t(i)-15)^4/(1+(t(i)-15)^4);
   
   dMsF = b*FSH(i) + (c1*FSH(i) - c2*LH(i)^a)*MsF(i-1);
   MsF(i) = MsF(i-1) + h * dMsF;
   dGrF = c2*LH(i)^a*MsF(i-1) + c3*LH(i)^a - c4*LH(i)*GrF(i-1);
   GrF(i) = GrF(i-1) + h * dGrF;
   dPrF = c4*LH(i)*GrF(i-1) - c5*LH(i)^a*PrF(i-1);
   PrF(i) = PrF(i-1) + h*dPrF;
   dSc1 = c5*LH(i)^a*PrF(i-1) - d1* Sc1(i-1);
   Sc1(i) = Sc1(i-1) + h*dSc1;
   dSc2 = d1*Sc1(i-1) - d2* Sc2(i-1);
   Sc2(i) = Sc2(i-1) + h*dSc2; 
   dLut1 = d2*Sc2(i-1) - k1*Lut1(i-1);
   Lut1(i) = Lut1(i-1) + h*dLut1;
   dLut2 = k1*Lut1(i-1) - k2*Lut2(i-1);
   Lut2(i) = Lut2(i-1) + h*dLut2;
   dLut3 = k2*Lut2(i-1) - k3*Lut3(i-1);
   Lut3(i) = Lut3(i-1) + h*dLut3;   
   dLut4 = k3*Lut3(i-1) - k4*Lut4(i-1);
   Lut4(i) = Lut4(i-1) + h*dLut4;
   
   
   E2(i) = e0 + e1*GrF(i) + e2*PrF(i) + e3*Lut4(i); %53 nmol/L
   P4(i) = p1*Lut3(i) + p2*Lut4(i); %1 nmol/L
   Ih(i) = h0 + h1 * PrF(i) + h2 * Lut2(i) + h3 * Lut3(i) + h4 * Lut4(i); %290 U/L   
end
figure(1);
subplot(3,1,1); hold on; grid minor
plot(t, E2, 'Linewidth', 2)
xlabel('days')
ylabel('E_2')
title('The amount of Estradiol (E_2) in blood')

subplot(3,1,2); hold on; grid minor
plot(t, P4, 'Linewidth', 2)
xlabel('days')
ylabel('P_4')
title('The amount of Progesterone (P_4) in blood')

subplot(3,1,3); hold on; grid minor
plot(t, Ih, 'Linewidth', 2)
xlabel('days')
ylabel('I_h')
title('The amount of Inhibin (Ih) in blood')

sgtitle('Ovarian Hormones')

figure(2);
subplot(2,1,1); hold on; grid minor
plot(t, FSH, 'Linewidth', 2)
xlabel('days')
ylabel('FSH')
title('The amount of FSH in blood')

subplot(2,1,2); hold on; grid minor
plot(t, LH, 'Linewidth', 2)
xlabel('days')
ylabel('LH')
title('The amount of LH in blood')


sgtitle('Gonadotropin Hormones')