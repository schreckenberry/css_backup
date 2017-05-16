kappaAppr = csvread('/home/max/gits/css-new/build/Data/kappaAppr.csv');
kappaConv = csvread('/home/max/gits/css-new/build/Data/kappaConv.csv');

sigma=kappaAppr(:,1);
kappaAppr = kappaAppr(:,2:end);
kappaConv = kappaConv(:,2:end);

diff = abs(kappaAppr-kappaConv);
figure;
surf(1:size(diff,2),sigma,diff,'Edgecolor','none');
title('difference of kappa (appr/conv)');
xlabel('idx');
ylabel('sigma');
zlabel('diff');
view(160, 45);


mseNorm = 100*sqrt(diff.^2)./min(max(max(kappaAppr)),max(max(kappaConv)));
disp(['The mean MSE for the Approximation of Kappa is ' num2str(mean(mean(mseNorm))) '%']);
disp(['The max. MSE for the Approximation of Kappa is ' num2str(max(max(mseNorm))) '%']);



%%
figure;
subplot(2,2,1);
plot(kappaConv(1,:)); hold on;
plot(kappaAppr(1,:),'--');
legend('conv','appr');
title(['kappa, sigma=' num2str(sigma(1))]);

subplot(2,2,2);
plot(kappaConv(2,:)); hold on;
plot(kappaAppr(2,:),'--');
legend('conv','appr');
title(['kappa, sigma=' num2str(sigma(2))]);

subplot(2,2,3);
plot(kappaConv(3,:)); hold on;
plot(kappaAppr(3,:),'--');
legend('conv','appr');
title(['kappa, sigma=' num2str(sigma(3))]);

subplot(2,2,4);
plot(kappaConv(6,:)); hold on;
plot(kappaAppr(6,:),'--');
legend('conv','appr');
title(['kappa, sigma=' num2str(sigma(6))]);