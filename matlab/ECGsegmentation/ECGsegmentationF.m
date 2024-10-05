%% UNIVERSIDADE FEDERAL DO MARANHÃO
%% MESTRADO DE ENGENHARIA DE ELETRICIDADE
%% JONATHAN ARAUJO QUEIROZ

%%        CÓDIGO : ECG Segmentation
%%        DATA: 11/05/2018  


function [B,P,QRS,T] = ECGsegmentationF(signal,Fs)
%%%%%%%%%%%%%%%%%%%% Pametros %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('ECG.mat');
%signal=[signal(5000:100000,:)];
%Fs = 81;
Fs=round(Fs*0.7);
theta=0.4;
lambda=0.6;
lambdaP=0.15;
thetaQRS=0.15;
lambdaQRS=0.15;

%%%%%%%%%%%%%%%%%%%%%%% Pre-processamento %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal=signal(round(size(signal,1)*0.01):end-(round(size(signal,1)*0.01))); % Retirar 1% do inicio e fim por erro de amostragem
signal=signal-mean(signal); % retirar média
signal=signal ./(sdt(signal)); % Normalizar
[signal, qrs_amplitude,qrs_index,delay] = detection_peaksR(signal, [0 60],  Fs, 0); % detecção de picos 

%%%%%%%%%%%%%%%%%%%%%%% Ciclo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(qrs_index(1)<theta*Fs) %Exclui picos para garantir o intervalo;
    qrs_index(1)=[];
    qrs_amplitude(1)=[];
end

if(qrs_index(end)>length(signal)-lambda*Fs) %Exclui picos para garantir o intervalo;
    qrs_index(end)=[];
    qrs_amplitude(end)=[];
end
% % % % % % % % Só lembrando que o primeiro if, é pz\ara os 
% % % % % % % % pontos antes do pico e o segundo para os pontos depois

B = zeros(length(qrs_amplitude),Fs+1);

%%%%%%%%%%%%%%%%%%%%%%% Batimentos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor i = 1:length(qrs_amplitude)
    % Batimentos
    B(i,:) = signal(qrs_index(i)-round(theta*Fs):qrs_index(i)+round(lambda*Fs));
    % Onda P
    P(i,:) = signal(qrs_index(i)-round(theta*Fs):(qrs_index(i)-round(theta*Fs))+round(lambdaP*Fs));   
    % Complexo QRS
    QRS(i,:) = signal(qrs_index(i)-round(thetaQRS*Fs):qrs_index(i)+round(lambdaQRS*Fs));   
    % Onda T
    T(i,:) = signal((qrs_index(i)+round(lambda*Fs))-round(0.3*Fs):qrs_index(i)+round(lambda*Fs));   
end 
B=B';
P=P';
QRS=QRS';
T=T';

%%%%%%%%%%%%%%%%%%%% figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nf=16;
figure
subplot(2,1,1)
set(gca,'fontsize',nf)
%hold on
plot(B)
xlabel('Samples', 'FontSize', nf) ;
ylabel('mV', 'FontSize', nf) ;
title('Healthy beats (50000)', 'FontSize', nf) ;

subplot(2,3,4)
set(gca,'fontsize',nf)
%hold on
plot(P)
xlabel('Samples', 'FontSize', nf) ;
ylabel('mV', 'FontSize', nf) ;
title('Healthy beats (50000)', 'FontSize', nf) ;

subplot(2,3,5)
set(gca,'fontsize',nf)
%hold on
plot(QRS)
xlabel('Samples', 'FontSize', nf) ;
ylabel('mV', 'FontSize', nf) ;
title('Healthy beats (50000)', 'FontSize', nf) ;

subplot(2,3,6)
set(gca,'fontsize',nf)
%hold on
plot(T)
xlabel('Samples', 'FontSize', nf) ;
ylabel('mV', 'FontSize', nf) ;
title('Healthy beats (50000)', 'FontSize', nf) ;

end

