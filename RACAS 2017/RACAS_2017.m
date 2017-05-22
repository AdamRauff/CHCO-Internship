% RACAS 2017 Data

% boxplots by WHO-FC
load('RACAS_DATA.mat');
subplot(2,2,1);
boxplot(KHVVCR, whoclass);
title('VVCR by WHO Functional Class','FontSize',20);
xlabel('WHO functional class','FontSize',18);
ylabel('VVCR','FontSize',18);
hold off

subplot(2,2,2);
boxplot(pvr, whoclass);
title('PVR by WHO Functional Class','FontSize',20);
xlabel('WHO functional class','FontSize',18);
ylabel('PVR','FontSize',18);
hold off

% Boxplots by death
subplot(2,2,3);
boxplot(KHVVCR, death);
title('VVCR by death','FontSize',20);
xlabel('death','FontSize',18);
ylabel('VVCR','FontSize',18);
hold off

subplot(2,2,4);
boxplot(pvr, death);
title('PVR by death','FontSize',20);
xlabel('death','FontSize',18);
ylabel('PVR','FontSize',18);
hold off