%% UNIVERSIDADE FEDERAL DO MARANHÃO
%% MESTRADO DE ENGENHARIA DE ELETRICIDADE
%% JONATHAN ARAUJO QUEIROZ


%%        CÓDIGO : ECG Segmentation
%%        DATA: 30/05/2018  

%clear all
close all
clc
%%%%%%%%%%% Relevant Component Analysis %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Normal
load('B.mat')
meann=mean(B);
stdn=std(B);
varn=var(B);
skewnessn=skewness(B);
kurtosisn=kurtosis(B);
T_n=[kurtosisn;meann;skewnessn];
T_n1=[(kurtosisn+wentropy(kurtosisn,'shannon'));meann+wentropy(meann,'shannon');skewnessn+wentropy(skewnessn,'shannon')];


%%%%%%%%%%%% AF
load('Baf.mat')
meanaf=mean(Baf);
stdaf=std(Baf);
varaf=var(Baf);
skewnessaf=skewness(Baf);
kurtosisaf=kurtosis(Baf);
T_n=[kurtosisaf;meanaf;skewnessaf];
T_af1=[(kurtosisaf+wentropy(kurtosisaf,'shannon'));meanaf+wentropy(meanaf,'shannon');skewnessaf+wentropy(skewnessaf,'shannon')];

%%%%%%%%%% Arritmias
load('Ba.mat')
meana=mean(Ba);
stda=std(Ba);
vara=var(Ba);
skewnessa=skewness(Ba);
kurtosisa=kurtosis(Ba);
T_a=[kurtosisa;meana;skewnessa];
T_a1=[(kurtosisa+wentropy(kurtosisa,'shannon'));meana+wentropy(meana,'shannon');skewnessa+wentropy(skewnessa,'shannon')];

%%%%%%%%%% Apneia 
load('Bpa1.mat')
meanap1=mean(Bap1);
stdap1=std(Bap1);
varap1=var(Bap1);
skewnessap1=skewness(Bap1);
kurtosisap1=kurtosis(Bap1);
T_ap1=[kurtosisap1;meanap1;skewnessap1];
T_ap1=[(kurtosisap1+wentropy(kurtosisap1,'shannon'));meanap1+wentropy(meanap1,'shannon');skewnessap1+wentropy(skewnessap1,'shannon')];

%%%%%%%%%% Apneia 2
load('Bap2.mat')
meanap2=mean(Bap2);
stdap2=std(Bap2);
varap2=var(Bap2);
skewnessap2=skewness(Bap2);
kurtosisap2=kurtosis(Bap2);
T_ap2=[kurtosisap2;meanap2;skewnessap2];
T_ap12=[(kurtosisap2+wentropy(kurtosisap2,'shannon'));meanap2+wentropy(meanap2,'shannon');skewnessap2+wentropy(skewnessap2,'shannon')];

%%%%%%%%%% Apneia 3
load('Bap3.mat')
meanap3=mean(Bap3);
stdap3=std(Bap3);
varap3=var(Bap3);
skewnessap3=skewness(Bap3);
kurtosisap3=kurtosis(Bap3);
T_ap3=[kurtosisap3;meanap3;skewnessap3];
T_ap13=[(kurtosisap3+wentropy(kurtosisap3,'shannon'));meanap3+wentropy(meanap3,'shannon');skewnessap3+wentropy(skewnessap3,'shannon')];


%%%%%%%%%% Supraventricular Arrhythmia
load('Bsa.mat')
meansa=mean(Bsa);
stdsa=std(Bsa);
varsa=var(Bsa);
skewnesssa=skewness(Bsa);
kurtosissa=kurtosis(Bsa);
T_sa=[kurtosissa;meansa;skewnesssa];
T_sa=[(kurtosissa+wentropy(kurtosissa,'shannon'));meansa+wentropy(meansa,'shannon');skewnesssa+wentropy(skewnesssa,'shannon')];


%%%%%%%%%% Morte Cardíaca
load('Bmc.mat')
meanmc=mean(Bmc);
stdmc=std(Bmc);
varmc=var(Bmc);
skewnessmc=skewness(Bmc);
kurtosismc=kurtosis(Bmc);
T_mc=[kurtosismc;meanmc;skewnessmc];
T_mc=[(kurtosismc+wentropy(kurtosismc,'shannon'));meanmc+wentropy(meanmc,'shannon');skewnessmc+wentropy(skewnessmc,'shannon')];



%%%%%%%%%% Insuficiência cardíaca congestiva
load('Bicc.mat')
meanicc=mean(Bicc);
stdicc=std(Bicc);
varicc=var(Bicc);
skewnessicc=skewness(Bicc);
kurtosisicc=kurtosis(Bicc);
T_icc=[kurtosisicc;meanicc;skewnessicc];
T_icc=[(kurtosisicc+wentropy(kurtosisicc,'shannon'));meanicc+wentropy(meanicc,'shannon');skewnessicc+wentropy(skewnessicc,'shannon')];




%%%%%%%%%% Epilepsy Post-Ictal
load('Bep.mat')
meanep=mean(Bep);
stdep=std(Bep);
varep=var(Bep);
skewnessep=skewness(Bep);
kurtosisep=kurtosis(Bep);
T_ep=[kurtosisep;meanep;skewnessep];
T_ep=[(kurtosisep+wentropy(kurtosisep,'shannon'));meanep+wentropy(meanep,'shannon');skewnessep+wentropy(skewnessep,'shannon')];



%%%%%%%%%% Abdominal and Direct Fetal
load('Badf.mat')
meanadf=mean(Badf);
stdadf=std(Badf);
varadf=var(Badf);
skewnessadf=skewness(Badf);
kurtosisadf=kurtosis(Badf);
T_adf=[kurtosisadf;meanadf;skewnessadf];
T_adf=[(kurtosisadf+wentropy(kurtosisadf,'shannon'));meanadf+wentropy(meanadf,'shannon');skewnessadf+wentropy(skewnessadf,'shannon')];




figure
subplot(2,2,1)
plot3(T_n1(1,:),T_n1(2,:),T_n1(3,:), 'y*', 'linewidth',2)
hold on
plot3(T_af1(1,:),T_af1(2,:),T_af1(3,:), 'm*', 'linewidth',2)
plot3(T_a1(1,:),T_a1(2,:),T_a1(3,:), 'c*', 'linewidth',2)
plot3(T_ap1(1,:),T_ap1(2,:),T_ap1(3,:), 'r*', 'linewidth',2)
plot3(T_ap12(1,:),T_ap12(2,:),T_ap12(3,:), 'g*', 'linewidth',2)
plot3(T_ap13(1,:),T_ap13(2,:),T_ap13(3,:), 'b*', 'linewidth',2)
plot3(T_sa(1,:),T_sa(2,:),T_sa(3,:), 'r*', 'linewidth',2)
plot3(T_mc(1,:),T_mc(2,:),T_mc(3,:), 'k*', 'linewidth',2)
plot3(T_icc(1,:),T_icc(2,:),T_icc(3,:), 'y+', 'linewidth',2)
plot3(T_ep(1,:),T_ep(2,:),T_ep(3,:), 'm+', 'linewidth',2)
plot3(T_adf(1,:),T_adf(2,:),T_adf(3,:), 'c+', 'linewidth',2)
xlabel({'Entropia+kurtosis'});
ylabel('Entropia+mean') ;
zlabel('Entropia+skewness');
legend('Healthy','AF','Arrhythmia','Apneia1','Apneia2','Apneia3','Supraventricular Arrhythmia','Morte Cardiaca','Insuficiencia cardiaca congestiva','Epilepsy Post-Ictal','Abdominal and Direct Fetal')



subplot(2,2,3)
plot3(T_n1(1,:),T_n1(2,:),T_n1(3,:), 'y*', 'linewidth',2)
hold on
plot3(T_af1(1,:),T_af1(2,:),T_af1(3,:), 'm*', 'linewidth',2)
plot3(T_a1(1,:),T_a1(2,:),T_a1(3,:), 'c*', 'linewidth',2)
plot3(T_ap1(1,:),T_ap1(2,:),T_ap1(3,:), 'r*', 'linewidth',2)
plot3(T_ap12(1,:),T_ap12(2,:),T_ap12(3,:), 'g*', 'linewidth',2)
plot3(T_ap13(1,:),T_ap13(2,:),T_ap13(3,:), 'b*', 'linewidth',2)
plot3(T_sa(1,:),T_sa(2,:),T_sa(3,:), 'r*', 'linewidth',2)
xlabel({'Entropia+kurtosis'});
ylabel('Entropia+mean') ;
zlabel('Entropia+skewness');
%legend('Healthy','AF','Arrhythmia','Apneia1','Apneia2','Apneia3','Supraventricular Arrhythmia')

subplot(2,2,4)
plot3(T_mc(1,:),T_mc(2,:),T_mc(3,:), 'k*', 'linewidth',2)
hold on
plot3(T_icc(1,:),T_icc(2,:),T_icc(3,:), 'y+', 'linewidth',2)
plot3(T_ep(1,:),T_ep(2,:),T_ep(3,:), 'm+', 'linewidth',2)
plot3(T_adf(1,:),T_adf(2,:),T_adf(3,:), 'c+', 'linewidth',2)
xlabel({'Entropia+kurtosis'});
ylabel('Entropia+mean') ;
zlabel('Entropia+skewness');
%legend('Morte Cardiaca','Insuficiencia cardiaca congestiva','Epilepsy Post-Ictal','Abdominal and Direct Fetal')
%title('593647 Heartbeats')




figure
plot3(T_ap1(1,:),T_ap1(2,:),T_ap1(3,:), 'r*', 'linewidth',2)
hold on
plot3(T_ap12(1,:),T_ap12(2,:),T_ap12(3,:), 'g*', 'linewidth',2)
plot3(T_ap13(1,:),T_ap13(2,:),T_ap13(3,:), 'b*', 'linewidth',2)


% figure
% plot3(T_n1(1,:),T_n1(2,:),T_n1(3,:), 'y*', 'linewidth',2)
% hold on
% plot3(T_af1(1,:),T_af1(2,:),T_af1(3,:), 'm*', 'linewidth',2)
% plot3(T_a1(1,:),T_a1(2,:),T_a1(3,:), 'c*', 'linewidth',2)
% plot3(T_ap1(1,:),T_ap1(2,:),T_ap1(3,:), 'r*', 'linewidth',2)
% plot3(T_ap12(1,:),T_ap12(2,:),T_ap12(3,:), 'g*', 'linewidth',2)
% plot3(T_ap13(1,:),T_ap13(2,:),T_ap13(3,:), 'b*', 'linewidth',2)
% plot3(T_sa(1,:),T_sa(2,:),T_sa(3,:), 'w*', 'linewidth',2)
% plot3(T_mc(1,:),T_mc(2,:),T_mc(3,:), 'k*', 'linewidth',2)
% plot3(T_icc(1,:),T_icc(2,:),T_icc(3,:), 'y+', 'linewidth',2)
% plot3(T_ep(1,:),T_ep(2,:),T_ep(3,:), 'm+', 'linewidth',2)
% plot3(T_adf(1,:),T_adf(2,:),T_adf(3,:), 'c+', 'linewidth',2)
% xlabel({'Entropia+kurtosis'});
% ylabel('Entropia+mean') ;
% zlabel('Entropia+skewness');
% legend('Healthy','AF','Arrhythmia','Apneia1','Apneia2','Apneia3','Supraventricular Arrhythmia','Morte Cardiaca','Insuficiencia cardiaca congestiva','Epilepsy Post-Ictal','Abdominal and Direct Fetal')
% title('593647 Heartbeats')

break



figure
subplot(2,4,1)
plot3(T_n(1,:),T_n(2,:),T_n(3,:), 'b*', 'linewidth',2)
xlabel({'Kurtosis'}) ;
ylabel('Mean') ;
zlabel('Skewness') ;
title('192222 Heartbeats')
%legend('Healthy')
%legend('Healthy','Arrhythmia','Apneia')

subplot(2,4,2)
plot3(T_a(1,:),T_a(2,:),T_a(3,:), 'r*', 'linewidth',2)
xlabel({'Kurtosis'}) ;
ylabel('Mean') ;
zlabel('Skewness') ;
title('175107 Heartbeats')
%legend('Arrhythmia')

subplot(2,4,3)
plot3(T_ap(1,:),T_ap(2,:),T_ap(3,:), 'm*', 'linewidth',2)
xlabel({'Kurtosis'}) ;
ylabel('Mean') ;
zlabel('Skewness') ;
title('226318 Heartbeats')

subplot(2,4,4)
plot3(T_n(1,:),T_n(2,:),T_n(3,:), 'b*', 'linewidth',2)
hold on
plot3(T_a(1,:),T_a(2,:),T_a(3,:), 'r*', 'linewidth',2)
plot3(T_ap(1,:),T_ap(2,:),T_ap(3,:), 'm*', 'linewidth',2)
xlabel({'Kurtosis'}) ;
ylabel('Mean') ;
zlabel('Skewness') ;
legend('Healthy','Arrhythmia','Apneia')
title('593647 Heartbeats')

%%%%%
subplot(2,4,5)
plot3(T_n1(1,:),T_n1(2,:),T_n1(3,:), 'b*', 'linewidth',2)
xlabel({'Entropia+kurtosis'});
ylabel('Entropia+mean') ;
zlabel('Entropia+skewness');
title('192222 Heartbeats')
%legend('Healthy')
%legend('Healthy','Arrhythmia','Apneia')

subplot(2,4,6)
plot3(T_a1(1,:),T_a1(2,:),T_a1(3,:), 'r*', 'linewidth',2)
xlabel({'Entropia+kurtosis'});
ylabel('Entropia+mean') ;
zlabel('Entropia+skewness');
title('175107 Heartbeats')
%legend('Arrhythmia')

subplot(2,4,7)
plot3(T_ap1(1,:),T_ap1(2,:),T_ap1(3,:), 'm*', 'linewidth',2)
xlabel({'Entropia+kurtosis'});
ylabel('Entropia+mean') ;
zlabel('Entropia+skewness');
title('226318 Heartbeats')

subplot(2,4,8)
plot3(T_n1(1,:),T_n1(2,:),T_n1(3,:), 'b*', 'linewidth',2)
hold on
plot3(T_af1(1,:),T_af1(2,:),T_af1(3,:), 'y*', 'linewidth',2)
plot3(T_a1(1,:),T_a1(2,:),T_a1(3,:), 'r*', 'linewidth',2)
plot3(T_ap1(1,:),T_ap1(2,:),T_ap1(3,:), 'm*', 'linewidth',2)
plot3(T_ap12(1,:),T_ap12(2,:),T_ap12(3,:), 'g*', 'linewidth',2)
plot3(T_ap13(1,:),T_ap13(2,:),T_ap13(3,:), 'k*', 'linewidth',2)
xlabel({'Entropia+kurtosis'});
ylabel('Entropia+mean') ;
zlabel('Entropia+skewness');
legend('Healthy','Arrhythmia','Apneia1','Apneia2','Apneia3')
title('593647 Heartbeats')
%grid on


figure
subplot(2,4,1)
histfit(kurtosisn,500,'kernel')
subplot(2,4,2)
histfit(kurtosisa,500,'kernel')
subplot(2,4,3)
histfit(kurtosisap,500,'kernel')
subplot(2,4,4)
histfit(kurtosisn,500,'kernel')
hold on
histfit(kurtosisa,500,'kernel')
histfit(kurtosisap,500,'kernel')

subplot(2,4,5)
histfit(kurtosisn+wentropy(kurtosisn,'shannon'),500,'kernel')
subplot(2,4,6)
histfit(kurtosisa+wentropy(kurtosisa,'shannon'),500,'kernel')
subplot(2,4,7)
histfit(kurtosisap+wentropy(kurtosisap,'shannon'),500,'kernel')
subplot(2,4,8)
histfit(kurtosisn+wentropy(kurtosisn,'shannon'),500,'kernel')
hold on
histfit(kurtosisa+wentropy(kurtosisa,'shannon'),500,'kernel')
histfit(kurtosisap+wentropy(kurtosisap,'shannon'),500,'kernel')

break

%%        CÓDIGO : ECG Segmentation
%%        DATA: 22/05/2018  


%%%%%%%%%%% Relevant Component Analysis %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Normal
meann=mean(B);
stdn=std(B);
varn=var(B);
skewnessn=skewness(B);
kurtosisn=kurtosis(B);

RCAn=[kurtosisn;meann;skewnessn];

RCAn_ms=meann*skewnessn';
RCAn_sk=skewnessn*kurtosisn';
RCAn_km=kurtosisn*meann';

An1=[RCAn_ms 0 0]*eye(3);

T_n=[RCAn_ms 0 0;0 RCAn_sk 0;0 0 RCAn_km]*RCAn;


%%%%%%%%%% Arritmias
meana=mean(BA);
stda=std(BA);
vara=var(BA);
skewnessa=skewness(BA);
kurtosisa=kurtosis(BA);

RCAa=[kurtosisa;meana;skewnessa];

RCAa_ms=meana*skewnessa';
RCAa_sk=skewnessa*kurtosisa';
RCAa_km=kurtosisa*meana';

T_a=[RCAa_ms 0 0;0 RCAa_sk 0;0 0 RCAa_km]*RCAa;


%%%%%%%%%% Apneia 
meanap=mean(BAp);
stdap=std(BAp);
varap=var(BAp);
skewnessap=skewness(BAp);
kurtosisap=kurtosis(BAp);

RCAap=[kurtosisap;meanap;skewnessap];

RCAap_ms=meanap*skewnessap';
RCAap_sk=skewnessap*kurtosisap';
RCAap_km=kurtosisap*meanap';

T_ap=[RCAap_ms 0 0;0 RCAap_sk 0;0 0 RCAap_km]*RCAap;



figure
subplot(1,2,1)
plot3(meann,skewnessn,kurtosisn, 'b*', 'linewidth',2)
hold on 
plot3(meana,skewnessa,kurtosisa, 'r*', 'linewidth',2)
plot3(meanap,skewnessap,kurtosisap, 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
%legend('Healthy','Arrhythmia','Apneia')


subplot(1,2,2)
plot3(T_n(1,:),T_n(2,:),T_n(3,:), 'b*', 'linewidth',2)
hold on
plot3(T_a(1,:),T_a(2,:),T_a(3,:), 'r*', 'linewidth',2)
plot3(T_ap(1,:),T_ap(2,:),T_ap(3,:), 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
%legend('Healthy','Arrhythmia','Apneia')


break
% 
% % %%% \\\\\\\\\\\\\\ análise fatorial\\\\\\\\\\\\\\\\\\\\\\\\\
% % %%%%%%%%%%%% Normal
% % Itemn=[meann;stdn;varn;skewnessn;kurtosisn];
% % RCAn=factoran(Itemn,1);
% % 
% % %%%%%%%%%% Arritmias
% % Itema=[meana;stda;vara;skewnessa;kurtosisa];
% % RCAa=factoran(Itema,1);
% % 
% % 
% % %%%%%%%%%% Apneia 
% % Itemap=[meanap stdap varap skewnessap kurtosisap];
% % RCAap=factoran(Itemap,1);
% 
% 
% 
% %%% \\\\\\\\\\\\\\  \\\\\\\\\\\\\\\\\\\\\\\\\
% %%%%%%%%%%%% Normal
% RCAn=[meann;skewnessn;kurtosisn];
% 
% %%%%%%%%%% Arritmias
% RCAa=[meana;skewnessa;kurtosisa];
% 
% 
% %%%%%%%%%% Apneia 
% Itemap=[meanap stdap varap skewnessap kurtosisap];
% RCAap=[meanap;skewnessap;kurtosisap];


%%%%%% fig artigo
x0 = ecg(64);
x1 = ecg(128);
x2 = ecg(256);
x3 = ecg(512);
x4 = ecg(1024);
x5 = ecg(2048);
x6 = ecg(4096);


figure

subplot(4,2,1)
plot(x1);
hold on; 
plot(x1,'r*');
subplot(4,4,3)
histfit(x1,10)
subplot(4,4,4)
boxplot(x1)

mean0=mean(x0);
std0=std(x0);
var0=var(x0);
skewness0=skewness(x0);
kurtosis0=kurtosis(x0);
w0=wentropy(x0,'shannon');

mean1=mean(x1);
std1=std(x1);
var1=var(x1);
skewness1=skewness(x1);
kurtosis1=kurtosis(x1);
w1=wentropy(x1,'shannon');

mean2=mean(x2);
std2=std(x2);
var2=var(x2);
skewness2=skewness(x2);
kurtosis2=kurtosis(x2);
w2=wentropy(x2,'shannon');

mean3=mean(x3);
std3=std(x3);
var3=var(x3);
skewness3=skewness(x3);
kurtosis3=kurtosis(x3);
w3=wentropy(x3,'shannon');

mean4=mean(x4);
std4=std(x4);
var4=var(x4);
skewness4=skewness(x4);
kurtosis4=kurtosis(x4);
w4=wentropy(x4,'shannon');


mean5=mean(x5);
std5=std(x5);
var5=var(x5);
skewness5=skewness(x5);
kurtosis5=kurtosis(x5);
w5=wentropy(x5,'shannon');


mean6=mean(x6);
std6=std(x6);
var6=var(x6);
skewness6=skewness(x6);
kurtosis6=kurtosis(x6);
w6=wentropy(x6,'shannon');


subplot(4,2,3)
plot(x2);
hold on; 
plot(x2,'b*');
subplot(4,4,7)
histfit(x2,10)
subplot(4,4,8)
boxplot(x2)


subplot(4,2,5)
plot(x3);
hold on; 
plot(x3,'g*');
subplot(4,4,11)
histfit(x3,10)
subplot(4,4,12)
boxplot(x3)

subplot(4,2,7)
plot(x4);
hold on; 
plot(x4,'m*');
subplot(4,4,15)
histfit(x4,10)
subplot(4,4,16)
boxplot(x4)


