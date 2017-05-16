% apprM = csvread('/home/max/gits/css-or/build/Data/approximation.csv');
% convM = csvread('/home/max/gits/css-or/build/Data/convolution.csv');

apprM = csvread('/home/max/gits/css-new/build/Data/appr0.csv');
convM = csvread('/home/max/gits/css-new/build/Data/conv0.csv');

% apprM = csvread('/home/max/gits/css-new/build/Data/appr1.csv');
% convM = csvread('/home/max/gits/css-new/build/Data/conv1.csv');

% apprM = csvread('/home/max/gits/css-new/build/Data/appr2.csv');
% convM = csvread('/home/max/gits/css-new/build/Data/conv2.csv');

sigma=8:2:200;

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

%% MSE Calculation
% I have no idea what's the best check is...
diff=abs(apprM-convM);
mse = mean(diff,2);
msepercent = 100*mse./max(max(convM));

figure;
plot(sigma, msepercent);
xlabel('sigma');
ylabel('MSE [%]');
title('MSE of Approximation for each sigma');

disp(['worst absolute MSE=' num2str(max(mse))]);
disp(['worst relative MSE=' num2str(max(msepercent)) '%']);
disp(['mean absolute MSE=' num2str(mean(mse))]);
disp(['mean relative MSE=' num2str(mean(msepercent)) '%']);


