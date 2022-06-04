clear;clc;close all
h = 0.01;
t = 0:h:30;

for i = 1:length(t)
   val = 175*exp(-(t(i)-5)^2/120) + 150*exp(-(t(i)-35)^2/160);
   FSH(i) = 250 - (250*(t(i) - 15)^4 )/(1+(t(i)-15)^4) + val;
   LH(i) = 380 - 352*(t(i)-15)^4/(1+(t(i)-15)^4);
end
figure(1);
subplot(2,1,1); hold on; grid minor
plot(t, FSH, 'Linewidth', 2)
xlabel('days')
ylabel('FSH [mg/L]')
title('The amount of FSH in blood')

subplot(2,1,2); hold on; grid minor
plot(t, LH, 'Linewidth', 2)
xlabel('days')
ylabel('LH [mg/L]')
title('The amount of LH in blood')

sgtitle('Hormones')