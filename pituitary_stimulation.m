clc;clear;clear all;close all
%%
h = 0.01;
t = 0:h:30;
%%
dE = 0.42; %days
dP = 2.9; %days
v_dis = 2.5; %L
dIh = 2; %days
%%
% Ih(1) = 406;
% P4(1) = 3.418;
% E2(1) = 240.15;
Ih(1) = 0;
P4(1) = 0;
E2(1) = 0;

%%
dRPLH(1) = synLH(0, dE, dP) - relLH(0, 0);
RPLH(1) = 100 + h*dRPLH(1);
dLH(1) = 1/v_dis*relLH(0, RPLH(1)) - clearLH(0);
LH(1) = 100 + h*dLH(1);

dRPFSH(1) = synFSH(0, dIh) - relFSH(0, 0);
RPFSH(1) = 100 + h*dRPFSH(1);
dFSH(1) = 1/v_dis*relFSH(0, RPFSH(1)) - clearFSH(0);
FSH(1) = 100 + h* dFSH(1);
%%
for i = 2:length(t)

    E2(i) = E2FUNC(t(i), 0);
    P4(i) = P4FUNC(t(i), 0);
    Ih(i) = IhFUNC(t(i), 0);

    %%
%     newfunc(i) = 1.4+(95.9*1e+6*E2FUNC(t(i), 0)^8)/(360^8+E2FUNC(t(i), 0)^8);
    %%
    dRPLH(i) = synLH(t(i), dE, dP) - relLH(t(i), RPLH(i-1));
    RPLH(i) = RPLH(i-1) + h*dRPLH(i);
    dLH(i) = 1/v_dis*relLH(t(i), RPLH(i)) - clearLH(LH(i-1));
    LH(i) = LH(i-1) + h*dLH(i);

    dRPFSH(i) = synFSH(t(i), dIh) - relFSH(t(i), RPFSH(i-1));
    RPFSH(i) = RPFSH(i-1) + h*dRPFSH(i);
    dFSH(i) = 1/v_dis*relFSH(t(i), RPFSH(i)) - clearFSH(FSH(i-1));
    FSH(i) = FSH(i-1) + h* dFSH(i);
   
end
%%

figure(1)
subplot(3,1,1)
hold on; grid minor;
plot(t, E2)
title("Estradiol (E_2)")

subplot(3,1,2)
hold on; grid minor;
plot(t, P4)
title('Progesterone (P_4)')

subplot(3,1,3)
plot(t, Ih)
hold on; grid minor;
title('Inhibin (Ih)')
xlabel('time (s)')
%%
figure(2)
subplot(2,2,1)
hold on; grid minor;
plot(t, RPLH)
title("RPLH")
xlabel('time (s)')

subplot(2,2,2)
hold on; grid minor;
plot(t, LH)
title('LH')
xlabel('time (s)')

subplot(2,2,3)
plot(t, RPFSH)
hold on; grid minor;
title('RPFSH')
xlabel('time (s)')

subplot(2,2,4)
plot(t, FSH)
hold on; grid minor;
title('FSH')

xlabel('time (s)')
sgtitle('Hormones')
%%
% figure(3)
% plot(E2, newfunc./1e+6)
