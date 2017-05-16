%% Script for checking the results of CSSApproximation::approximateConvGaussAll
%% Part 3 - second derivative of gaussian

% Approximation of Convolving with the second derivative of a gaussian
appr2ndM = csvread('/home/max/gits/css-new/build/Data/appr2.csv');
conv2ndM = csvread('/home/max/gits/css-new/build/Data/conv2.csv');
sigma=8:2:200;

diff2nd = abs(conv2ndM-appr2ndM);
mse2nd = sqrt(diff2nd.^2);
mseNorm2nd = 100*mse2nd./min(max(max(appr2ndM)),max(max(conv2ndM)));

figure;
subplot(2,2,1);
plot(conv2ndM(1,:)); hold on;
plot(appr2ndM(1,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(1))]);

subplot(2,2,2);
plot(conv2ndM(2,:)); hold on;
plot(appr2ndM(2,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(2))]);

subplot(2,2,3);
plot(conv2ndM(3,:)); hold on;
plot(appr2ndM(3,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(3))]);

subplot(2,2,4);
plot(conv2ndM(6,:)); hold on;
plot(appr2ndM(6,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(6))]);

disp(['The mean MSE for the Approximation of a Conv. w/ the 2nd derivative of a gaussian is ' num2str(mean(mean(mseNorm2nd))) '%']);
disp(['The max. MSE for the Approximation of a Conv. w/ the 2nd derivative of a gaussian is ' num2str(max(max(mseNorm2nd))) '%']);
