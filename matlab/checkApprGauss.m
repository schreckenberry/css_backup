%% Script for checking the results of CSSApproximation::approximateConvGaussAll
%% Part 1 - Gaussian

% Reading the results
% Approximation of Convolving with a gaussian
apprM = csvread('/home/max/gits/css-new/build/Data/appr0.csv');
convM = csvread('/home/max/gits/css-new/build/Data/conv0.csv');
sigma=8:2:200;

diff = abs(convM-apprM);
mse = sqrt(diff.^2);
mseNorm = 100*mse./min(max(max(apprM)),max(max(convM)));

figure;
subplot(2,2,1);
plot(convM(10,:)); hold on;
plot(apprM(10,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(10))]);

subplot(2,2,2);
plot(convM(30,:)); hold on;
plot(apprM(30,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(30))]);

subplot(2,2,3);
plot(convM(60,:)); hold on;
plot(apprM(60,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(60))]);

subplot(2,2,4);
plot(convM(90,:)); hold on;
plot(apprM(90,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(90))]);

disp(['The mean MSE for the Approximation of a Conv. w/ a gaussian is ' num2str(mean(mean(mseNorm))) '%']);
disp(['The max. MSE for the Approximation of a Conv. w/ a gaussian is ' num2str(max(max(mseNorm))) '%']);