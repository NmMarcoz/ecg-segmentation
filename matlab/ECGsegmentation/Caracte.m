%%%%%%%%%%% relevant component analysis %%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%% morfologia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Normal
meann=mean(B);
stdn=std(B);
varn=var(B);
skewnessn=skewness(B);
kurtosisn=kurtosis(B);


% Hn=[];
% Hs=[];
% Hl=[];
% 
% for i=1:size(B,2)
% soma=sum(B(:,i).^2);
% pdf=((B(:,i).^2)./soma).^0.5;
% B(:,i)=pdf;
% Hs=[Hs wentropy(B(:,i),'shannon')];
% Hn=[Hn wentropy(B(:,i),'norm',1.1)];
% Hl=[Hl wentropy(B(:,i),'log energy')];
% end


%%%%%%%%%%%% Arritmias
meana=mean(BA);
stda=std(BA);
vara=var(BA);
skewnessa=skewness(BA);
kurtosisa=kurtosis(BA);

% Hsa=[];
% Hna=[];
% Hla=[];
% for i=1:size(BA,2)
% soma=sum(BA(:,i).^2);
% pdf=((BA(:,i).^2)./soma).^0.5;
% BA(:,i)=pdf;
% Hsa=[Hsa wentropy(BA(:,i),'shannon')];
% Hna=[Hna wentropy(BA(:,i),'norm',1.1)];
% Hla=[Hla wentropy(BA(:,i),'log energy')];
% end


%%%%%%%%%%%%%%%%% Apneia 


meanap=mean(BAp);
stdap=std(BAp);
varap=var(BAp);
skewnessap=skewness(BAp);
kurtosisap=kurtosis(BAp);
%wentropyap=wentropy(BAp,'shannon');



t1=(meann*skewnessn')*kurtosisn;
t2=(skewnessn*kurtosisn')*meann;
t3=(kurtosisn*meann')*skewnessn;

t1a=(meana*skewnessa')*kurtosisa;
t2a=(skewnessa*kurtosisa')*meana;
t3a=(kurtosisa*meana')*skewnessa;

t1ap=(meanap*skewnessap')*kurtosisap;
t2ap=(skewnessap*kurtosisap')*meanap;
t3ap=(kurtosisap*meanap')*skewnessap;

figure
plot3(t1,t2,t3, 'b*', 'linewidth',2)
hold on 
plot3(t1a,t2a,t3a, 'r*', 'linewidth',2)
plot3(t1ap,t2ap,t3ap, 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
legend('Healthy','Arrhythmia','Apneia')


%%%%%%%%%%%%% PCA %%%%%%%%%%%%%%%%%%%%
Xn=[meann;skewnessn;kurtosisn];
Xa=[meana;skewnessa;kurtosisa];
Xap=[meanap;skewnessap;kurtosisap];

[COEFFn,SCOREn] = princomp(Xn');
[COEFFa,SCOREa] = princomp(Xa');
[COEFFap,SCOREap] = princomp(Xap');



Txn=Xn'*COEFFn;
Txa=Xa'*COEFFa;
Txap=Xap'*COEFFap;

figure
plot3(meann,skewnessn,kurtosisn, 'b*', 'linewidth',2)
hold on 
plot3(meana,skewnessa,kurtosisa, 'r*', 'linewidth',2)
plot3(meanap,skewnessap,kurtosisap, 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
legend('Healthy','Arrhythmia','Apneia')



figure
plot3(Txn(:,1),Txn(:,2),Txn(:,3), 'b*', 'linewidth',2)
hold on 
plot3(Txa(:,1),Txa(:,2),Txa(:,3), 'r*', 'linewidth',2)
plot3(Txap(:,1),Txap(:,2),Txap(:,3), 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
legend('Healthy','Arrhythmia','Apneia')


break

% Hnp=[];
% Hsp=[];
% Hlp=[];
% 
% for i=1:size(BAp,2)
% soma=sum(BAp(:,i).^2);
% pdf=((BAp(:,i).^2)./soma).^0.5;
% BAp(:,i)=pdf;
% Hsp=[Hsp wentropy(BAp(:,i),'shannon')];
% Hnp=[Hnp wentropy(BAp(:,i),'norm',1.1)];
% Hlp=[Hlp wentropy(BAp(:,i),'log energy')];
% end

Nsts=stdn*skewnessn';
Nstk=stdn*kurtosisn';
Nsk=kurtosisn*skewnessn';

at=[Nsk;Nstk;Nsts];

Tn=[stdn;skewnessn;kurtosisn];
An=[Nsk Nstk Nsts; Nsk Nstk Nsts; Nsk Nstk Nsts];
Tnn=An'*Tn;


Asts=stda*skewnessa';
Astk=stda*kurtosisa';
Ask=kurtosisa*skewnessa';

Ta=[stda;skewnessa;kurtosisa];
Aa=[Ask Astk Asts; Ask Astk Asts; Ask Astk Asts];
Taa=Aa'*Ta;



Apsts=stdap*skewnessap';
Apstk=stdap*kurtosisap';
Apsk=kurtosisap*skewnessap';

Tap=[stdap;skewnessap;kurtosisap];
Aap=[Apsk Apstk Apsts; Apsk Apstk Apsts; Apsk Apstk Apsts];
Tapp=Aap'*Tap;


figure
plot3(Tnn(1,:),Tnn(2,:),Tnn(3,:), 'b*', 'linewidth',2)
hold on 
plot3(Taa(1,:),Taa(2,:),Taa(3,:), 'r*', 'linewidth',2)
plot3(Tapp(1,:),Tapp(1,:),Tapp(1,:), 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
legend('Healthy','Arrhythmia','Apneia')



break


Tn=[kurtosisn;skewnessn;stdn];
An=[(stdn*skewnessn') (stdn*kurtosisn') (kurtosisn*skewnessn');
    (stdn*skewnessn') (stdn*kurtosisn') (kurtosisn*skewnessn');
    (stdn*skewnessn') (stdn*kurtosisn') (kurtosisn*skewnessn')];
%An=An';
TAn=An'*Tn;

TA=[kurtosisa;skewnessa;stda];
AA=[(stda*skewnessa') (stda*kurtosisa') (kurtosisa*skewnessa');
    (stda*skewnessa') (stda*kurtosisa') (kurtosisa*skewnessa');
    (stda*skewnessa') (stda*kurtosisa') (kurtosisa*skewnessa')];
TAA=AA'*TA;


TAp=[kurtosisap;skewnessap;stdap];
AAp=[(stdap*skewnessap') (stdap*kurtosisap') (kurtosisap*skewnessap');
    (stdap*skewnessap') (stdap*kurtosisap') (kurtosisap*skewnessap');
    (stdap*skewnessap') (stdap*kurtosisap') (kurtosisap*skewnessap')];
TAAp=AAp'*TAp;


figure
plot3(TAn(1,:),TAn(2,:),TAn(3,:), 'b*', 'linewidth',2)
hold on 
plot3(TAA(1,:),TAA(2,:),TAA(3,:), 'r*', 'linewidth',2)
plot3(TAAp(1,:),TAAp(1,:),TAAp(1,:), 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
legend('Healthy','Arrhythmia','Apneia')


break
t1=(meann*skewnessn')*kurtosisn;
t2=(skewnessn*kurtosisn')*meann;
t3=(kurtosisn*meann')*skewnessn;

t1a=(meana*skewnessa')*kurtosisa;
t2a=(skewnessa*kurtosisa')*meana;
t3a=(kurtosisa*meana')*skewnessa;

t1ap=(meanap*skewnessap')*kurtosisap;
t2ap=(skewnessap*kurtosisap')*meanap;
t3ap=(kurtosisap*meanap')*skewnessap;

figure
plot3(t1,t2,t3, 'b*', 'linewidth',2)
hold on 
plot3(t1a,t2a,t3a, 'r*', 'linewidth',2)
plot3(t1ap,t2ap,t3ap, 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
legend('Healthy','Arrhythmia','Apneia')

x1=t1(:,1:size(t1a,2));
x2=t1a(:,1:size(t1a,2));
x3=t1ap(:,1:size(t1a,2));
x1=x1';
x2=x2';
x3=x3';
figure
boxplot([x1,x2,x3],{'Healthy','Arrhythmia','Apneia'})

figure
subplot(1,2,1)
plot3(t1,t2,t3, 'b*', 'linewidth',2)
hold on 
plot3(t1a,t2a,t3a, 'r*', 'linewidth',2)
plot3(t1ap,t2ap,t3ap, 'm*', 'linewidth',2)
xlabel({'MSK'}) ;
ylabel('SKM') ;
zlabel('KMS') ;
legend('Healthy','Arrhythmia','Apneia')

subplot(1,3,3)
boxplot([x1,x2,x3],{'Healthy','Arrhythmia','Apneia'})






break


figure
plot(Hs,'b*')
hold on
plot(Hsa,'r*')

figure
plot(Hn,'b*')
hold on
plot(Hna,'r*')

figure
plot(Hl,'b*')
hold on
plot(Hla,'r*')

%%%%%%%%%% nuvem


t1=meann*skewnessn'*kurtosisn;
t2=skewnessn*kurtosisn'*meann;
t3=kurtosisn*meann'*skewnessn;

t1a=meana*skewnessa'*kurtosisa;
t2a=skewnessa*kurtosisa'*meana;
t3a=kurtosisa*meana'*skewnessa;


figure
plot3(t1,t2,t3, 'b*', 'linewidth',2)
hold on 
plot3(t1a,t2a,t3a, 'r*', 'linewidth',2)



figure
plot3(meann,skewnessn,kurtosisn, 'b*', 'linewidth',2)
hold on 
plot3(meann,skewnessa,kurtosisa, 'r*', 'linewidth',2)


figure 
plot(Ns,'b+')
hold on
plot(Nsa,'m*')
xlabel('a) Samples') ;
ylabel('Average') ;
legend('Arrhythmia','Healthy');

figure 
plot(meann,'b+')
hold on
plot(meana,'m*')
xlabel('a) Samples') ;
ylabel('Average') ;
legend('Arrhythmia','Healthy');


figure
plot(vara,'m+')
hold on
plot(varn,'b*')
xlabel('a) Samples') ;
ylabel('Variance') ;
legend('Arrhythmia','Healthy');

figure
plot(skewnessa,'m+')
hold on
plot(skewnessn,'b*')
xlabel('a) Samples') ;
ylabel('Skewness') ;
legend('Arrhythmia','Healthy');


figure
plot(kurtosisa,'m+')
hold on
plot(kurtosisn,'b*')
xlabel('a) Samples') ;
ylabel('kurtosis') ;
legend('Arrhythmia','Healthy');



figure
plot3(mean,skewnessn,kurtosisn, 'b*', 'linewidth',2)
hold on 
plot3(meana,skewnessa,kurtosisa, 'r*', 'linewidth',2)
plot3(meanap,skewnessap,kurtosisap, 'y*', 'linewidth',2)
plot3(meanap1,skewnessap1,kurtosisap1, 'm*', 'linewidth',2)
xlabel({'Variance','b) Beats (100000)'}) ;
ylabel('Skewness') ;
zlabel('kurtosis') ;
legend('Healthy','Arrhythmia')


break

stdn=std(B);
varn=var(B);
skewnessn=skewness(B);
kurtosisn=kurtosis(B);

Item1=[meann;skewnessn;kurtosisn];
 
Anilfator1=factoran(Item1,1);

Item2=[meana stda vara skewnessa kurtosisa];
 
Anilfator2=factoran(Item2,1);



%  
% 
% %%%%%%%%%%%%%%%%%%%%%%% morfologia %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Normal
% meann=mean(B);
% stdn=std(B);
% varn=var(B);
% skewnessn=skewness(B);
% kurtosisn=kurtosis(B);

Hn=[];
Hs=[];
Hl=[];

for i=1:size(B,2)
% soma=sum(B(:,i).^2);
% pdf=((B(:,i).^2)./soma).^0.5;
% B(:,i)=pdf;
% Hs=[Hs wentropy(B(:,i),'shannon')];
% Hn=[Hn wentropy(B(:,i),'norm',1.1)];
% Hl=[Hl wentropy(B(:,i),'log energy')];

%%%%%%%%%%%%%%%% freu

% p = hist(BA);% calculate histogram counts
% p(p==0) = [];% remove zero entries in p
% BA = (p ./ length(BA)).^0.5;% normalize p so that sum(p) is one.
Hsf=[Hsf wentropy(B(:,i),'shannon')];
Hnf=[Hnf wentropy(B(:,i),'norm',1.1)];
Hlf=[Hnf wentropy(B(:,i),'log energy')];
end
                                 
                                 
%wentropyn=wentropy(B,'shannon');


% %break
% SIGMAn=[];
% ordemont=4;

% for i=1:ordemont
%  SIGMAn = [SIGMAn;moment(btn,i)];   
% end


% Arritmias
%load('BA.mat')

% meana=mean(BA);
% stda=std(BA);
% vara=var(BA);
% skewnessa=skewness(BA);
% kurtosisa=kurtosis(BA);

Has=[];
Han=[];
Hal=[];
Hasf=[];
Hanf=[];
Half=[];

for i=1:size(BA,2)
% soma=sum(BA(:,i).^2);
% pdf=((BA(:,i).^2)./soma).^0.5;
% BA(:,i)=pdf;
% Has=[Has wentropy(BA(:,i),'shannon')];
% Han=[Han wentropy(BA(:,i),'norm',1.1)];
% Hal=[Hal wentropy(BA(:,i),'log energy')];

%%%%%%%%%%%%%%%% freu

%p = hist(BA);% calculate histogram counts
%p(p==0) = [];% remove zero entries in p
%BA = (p ./ length(BA)).^0.5;% normalize p so that sum(p) is one.
% Hasf=[Hasf wentropy(BA(:,i),'shannon')];
% Hanf=[Hanf wentropy(BA(:,i),'norm',1.1)];
% Half=[Hanf wentropy(BA(:,i),'log energy')];
end

figure
plot(Hs,'b*')
hold on
plot(Has,'r*')

figure
plot(Hn,'b*')
hold on
plot(Han,'r*')

figure
plot(Hl,'b*')
hold on
plot(Hal,'r*')

plot3(Hs,Hn,Hl,'b*')
hold on
plot3(Has,Han,Hal,'r*')
break

% for i=1:size(BA,2)
% 
%  end                                

wentropya=wentropy(BA,'shannon');

meanap1=mean(BAp1);
stdap1=std(BAp1);
varap1=var(BAp1);
skewnessap1=skewness(BAp1);
kurtosisap1=kurtosis(BAp1);
wentropyap1=wentropy(BAp1,'shannon');

meanap=mean(BAp);
stdap=std(BAp);
varap=var(BAp);
skewnessap=skewness(BAp);
kurtosisap=kurtosis(BAp);
wentropyap=wentropy(BAp,'shannon');


%break
SIGMAa=[];
ordemonta=4;

% for i=1:ordemonta
%  SIGMAa = [SIGMAa;moment(bta,i)];   
% end




% %%%%%%%%%% nuvem
% figure 
% plot(meana,'m+')
% hold on
% plot(meann,'b*')
% xlabel('a) Samples') ;
% ylabel('Average') ;
% legend('Arrhythmia','Healthy');
% 
% 
% figure
% plot(vara,'m+')
% hold on
% plot(varn,'b*')
% xlabel('a) Samples') ;
% ylabel('Variance') ;
% legend('Arrhythmia','Healthy');
% 
% figure
% plot(skewnessa,'m+')
% hold on
% plot(skewnessn,'b*')
% xlabel('a) Samples') ;
% ylabel('Skewness') ;
% legend('Arrhythmia','Healthy');
% 
% 
% figure
% plot(kurtosisa,'m+')
% hold on
% plot(kurtosisn,'b*')
% xlabel('a) Samples') ;
% ylabel('kurtosis') ;
% legend('Arrhythmia','Healthy');

figure
plot3(meann,skewnessn,kurtosisn, 'b*', 'linewidth',2)
hold on 
plot3(meana,skewnessa,kurtosisa, 'r*', 'linewidth',2)
plot3(meanap,skewnessap,kurtosisap, 'y*', 'linewidth',2)
plot3(meanap1,skewnessap1,kurtosisap1, 'm*', 'linewidth',2)
xlabel({'Variance','b) Beats (100000)'}) ;
ylabel('Skewness') ;
zlabel('kurtosis') ;
legend('Healthy','Arrhythmia')