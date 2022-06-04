clear;clc;
h = 0.01;
t = 0:h:10000;
%%
% a = 0.6; b = 0.004; c1 = 0.0045; c2 = 0.077; c3 = 0.006; c4 = 0.008;
% c5 = 0.045; d1 = 0.5; d2 = 0.8; k1 = 0.6; k2 = 0.5; k3 = 0.6; k4 = 0.65;
%%
% e0 = 48; e1 = 0.7; e2 = 2.1; e3 = 1.7; p1 = 0.55; p2 = 0.45; h0 = 270;
% h1 = 2.5; h2 = 2; h3 = 10; h4 = 14;
%%
tau = 1.5; kLH = 2.42; alphaLH = 14; V0LH = 500; V1LH = 4500; KmLH = 200;
KiLHP = 12.2; CLHE = 0.004; CLHP = 0.26; VFSH = 375; alphaFSH = 8.21;
kFSH = 1.9; CFSHE = 0.0018; KiFSHIh = 3.5; CFSHP = 12; v = 2.5;
%%
b = 0.05; c1 = 0.08; c2 = 0.091; c3 = 0.13; c4 = 0.027; c5 = 0.51;
d1 = 0.5; d2 = 0.56; k1 = 0.55; k2 = 0.69; k3 = 0.85; k4 = 0.85;
alpha = 0.79; beta = 0.16; gama = 0.02; a = alpha;
%%
e0 = 30; e1 = 0.11; e2 = 0.21; e3 = 0.45; p0 = 0; p1 = 0.048; p2 = 0.048; h0 = 0.4;
h1 = 0.009; h2 = 0.029; h3 = 0.018;%; h4 = 14;

%%
   RPLH(1) = 29.65;
   LH(1) = 6.86;
   RPFSH(1) = 8.47;
   FSH(1) = 6.15;
   ReF(1) = 3.83;
   SeF(1) = 11.51;
   PrF(1) = 5.48;
   Ov1(1) = 19.27;
   Ov2(1) = 45.64;
   Lut1(1) = 100.73;
   Lut2(1) = 125.95;
   Lut3(1) = 135.84;   
   Lut4(1) = 168.71;

%    RPLH(1) = 40;
%    LH(1) = 12;
%    RPFSH(1) = 20;
%    FSH(1) = 11;
%    ReF(1) = 5;
%    SeF(1) = 1;
%    PrF(1) = 1;
%    Ov1(1) = 1;
%    Ov2(1) = 1;
%    Lut1(1) = 1;
%    Lut2(1) = 1;
%    Lut3(1) = 1;   
%    Lut4(1) = 1;
   
   E2(1) = e0 + e1*SeF(1) + e2*PrF(1) + e3*Lut4(1);
   P4(1) = p0 + p1*Lut3(1) + p2*Lut4(1);
   Ih(1) = h0 + h1 * PrF(1) + h2 * Lut3(1) + h3 * Lut4(1);
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
%    val = 175*exp(-(t(i)-5)^2/120) + 150*exp(-(t(i)-35)^2/160);
%    FSH(i) = 250 - (250*(t(i) - 15)^4 )/(1+(t(i)-15)^4) + val;
%    LH(i) = 380 - 352*(t(i)-15)^4/(1+(t(i)-15)^4);
   
%    Ih(i) = h0 + h1 * PrF(i) + h2 * Lut2(i) + h3 * Lut3(i) + h4 * Lut4(i); %290 U/L   
   
   dRPLHp1 = (V0LH+((V1LH*E2(i-1)^8)/(KmLH^8+E2(i-1)^8)))/(1+(P4(i-1))/KiLHP);
   dRPLHp2 = -(kLH*(1+CLHP*P4(i-1))*RPLH(i-1))/(1+CLHE*E2(i-1));
   
   dRPLH = dRPLHp1 + dRPLHp2;
   RPLH(i) = RPLH(i-1) + h * dRPLH;
   
   dLH = (1/v)*(kLH*(1+CLHP*P4(i-1))*RPLH(i-1))/(1+CLHE*E2(i-1))-alphaLH*LH(i-1);
   LH(i) = LH(i-1) + h * dLH;
   
   dRPFSH = VFSH/(1+(Ih(i-1)/KiFSHIh))-(kFSH*(1+CFSHP*P4(i-1))*RPFSH(i-1))/(1+CFSHE*E2(i-1)^2);
   RPFSH(i) = RPFSH(i-1) + h * dRPFSH;
   
   dFSH = (1/v)*((kFSH*(1+CFSHP*P4(i-1)))*RPFSH(i-1))/(1+CFSHE*E2(i-1)^2)-alphaFSH*FSH(i-1);
   FSH(i) = FSH(i-1) + h * dFSH;
   
%    dMsF = b*FSH(i) + (c1*FSH(i) - c2*LH(i)^a)*MsF(i-1);
%    MsF(i) = MsF(i-1) + h * dMsF;
   dReF = b*FSH(i-1) + (c1*FSH(i-1) - c2*LH(i-1)^a)*ReF(i-1);
   ReF(i) = ReF(i-1) + h * dReF;

   
   dSeF = c2*LH(i-1)^a*ReF(i-1) + (c3*LH(i-1)^beta - c4*LH(i-1))*SeF(i-1);
   SeF(i) = SeF(i-1) + h * dSeF;
   
   
   dPrF = c4*LH(i-1)*SeF(i-1) - c5*LH(i-1)^gama*PrF(i-1);
   PrF(i) = PrF(i-1) + h*dPrF;
   
   
%    dSc1 = c5*LH(i)^a*PrF(i-1) - d1* Sc1(i-1);
%    Sc1(i) = Sc1(i-1) + h*dSc1;
%    dSc2 = d1*Sc1(i-1) - d2* Sc2(i-1);
%    Sc2(i) = Sc2(i-1) + h*dSc2; 
   dOv1 = c5*LH(i-1)^gama*PrF(i-1) - d1* Ov1(i-1);
   Ov1(i) = Ov1(i-1) + h * dOv1;
   dOv2 = d1*Ov1(i-1) - d2* Ov2(i-1);
   Ov2(i) = Ov2(i-1) + h * dOv2;
   
   dLut1 = d2*Ov2(i-1) - k1*Lut1(i-1);
   Lut1(i) = Lut1(i-1) + h*dLut1;
   
   dLut2 = k1*Lut1(i-1) - k2*Lut2(i-1);
   Lut2(i) = Lut2(i-1) + h*dLut2;
   
   dLut3 = k2*Lut2(i-1) - k3*Lut3(i-1);
   Lut3(i) = Lut3(i-1) + h*dLut3;  
   
   dLut4 = k3*Lut3(i-1) - k4*Lut4(i-1);
   Lut4(i) = Lut4(i-1) + h*dLut4;
   
   E2(i) = e0 + e1*SeF(i-1) + e2*PrF(i-1) + e3*Lut4(i-1);
   P4(i) = p0 + p1*Lut3(i-1) + p2*Lut4(i-1);
   Ih(i) = h0 + h1 * PrF(i-1) + h2 * Lut3(i-1) + h3 * Lut4(i-1);
   
   y(1)=(V1LH*E2(1)^8)/(KmLH^8+E2(1)^8);
   y(i)=(V1LH*E2(i-1)^8)/(KmLH^8+E2(i-1)^8);
end
figure(1);
subplot(3,1,1); hold on; grid minor
plot(t, E2, 'Linewidth', 2)
xlabel('days')
ylabel('E_2 [ng/L]')
title('The amount of Estradiol (E_2) in blood')

subplot(3,1,2); hold on; grid minor
plot(t, P4, 'Linewidth', 2)
xlabel('days')
ylabel('P_4 [nmol/L]')
title('The amount of Progesterone (P_4) in blood')

subplot(3,1,3); hold on; grid minor
plot(t, Ih, 'Linewidth', 2)
xlabel('days')
ylabel('I_h [U/L]')
title('The amount of Inhibin (Ih) in blood')


sgtitle('Ovarian Hormones')
% figure(2);
% subplot(3,1,1); hold on; grid minor
% plot(t, ReF, 'Linewidth', 2)
% xlabel('days')
% ylabel('ReF')
% title('The amount of Estradiol (ReF) in blood')
% 
% subplot(3,1,2); hold on; grid minor
% plot(t, SeF, 'Linewidth', 2)
% xlabel('days')
% ylabel('SeF')
% title('The amount of Progesterone (SeF) in blood')
% 
% subplot(3,1,3); hold on; grid minor
% plot(t, PrF, 'Linewidth', 2)
% xlabel('days')
% ylabel('PrF')
% title('The amount of Inhibin (PrF) in blood')

figure(3);
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
% Km_LH = 50:h:400;
%%
% figure(4);hold on;grid on;
% plot(E2,y)
% 
% figure(5);hold on;grid on;
% plot(E2,LH)
