%% Script for checking the results of CSSApproximation::approximateConvGaussAll
%% Part 2 - first derivative of gaussian

% Approximation of Convolving with the first derivative of a gaussian
appr1stM = csvread('/home/max/gits/css-new/build/Data/appr1.csv');
conv1stM = csvread('/home/max/gits/css-new/build/Data/conv1.csv');
sigma=8:2:200;


diff1st = abs(conv1stM-appr1stM);
mse1st = sqrt(diff1st.^2);
mseNorm1st = 100*mse1st./min(max(max(appr1stM)),max(max(conv1stM)));

figure;
subplot(2,2,1);
plot(conv1stM(1,:)); hold on;
plot(appr1stM(1,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(1))]);

subplot(2,2,2);
plot(conv1stM(2,:)); hold on;
plot(appr1stM(2,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(2))]);

subplot(2,2,3);
plot(conv1stM(3,:)); hold on;
plot(appr1stM(3,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(3))]);

subplot(2,2,4);
plot(conv1stM(40,:)); hold on;
plot(appr1stM(40,:),'--');
legend('conv','appr');
title(['sigma=' num2str(sigma(40))]);

disp(['The mean MSE for the Approximation of a Conv. w/ the 1st derivative of a gaussian is ' num2str(mean(mean(mseNorm1st))) '%']);
disp(['The max. MSE for the Approximation of a Conv. w/ the 1st derivative of a gaussian is ' num2str(max(max(mseNorm1st))) '%']);