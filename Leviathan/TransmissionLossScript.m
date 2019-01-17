%This is a rough cut at including depth in the TL calculation
%In shallower water, the TL will more closely approximate the cylindrical
%curve

F = 10.^(-2:.1:1); %kHz
R = 1:10:40000; %yards
D1 = 4000; %deep - yards
D2 = 100; %shallow - yards

alpha = (.1*F.^2)./(1+F.^2) + (40*F.^2)./(4100+F.^2) + 2.75/10000 + 0.003;
logRmat = ones(size(F'))*log10(R);
alphaR = alpha' * R/1000;
TLs = 20*logRmat + alphaR;
TLc = 10*logRmat + alphaR;

% deep case
Rclose = 1:10:(D1/2); %range that is less than depth/2 (4000 yds - deep ocean)
Close = length(Rclose);
logRmat2 = ones(size(F'))*log10(Rclose);
TLsc1 = TLc; %start with calculated cylindrical spreading
TLsc1(:,1:Close) = TLsc1(:,1:Close)+ 10*logRmat2; %add 10logR for ranges less than depth
% add 10 log depth to account for the TL up to range = depth
TLsc1(:,(Close+1):end) = TLsc1(:,(Close+1):end)+ repmat(10*logRmat2(:,Close),1,(length(R)-Close)); 

% shallow case
Rclose = 1:10:(D2/2); %range that is less than depth/2 (100 yds - shallow ocean)
Close = length(Rclose);
logRmat2 = ones(size(F'))*log10(Rclose);
TLsc2 = TLc; %start with calculated cylindrical spreading
TLsc2(:,1:Close) = TLsc2(:,1:Close)+ 10*logRmat2; %add 10logR for ranges less than depth
% add 10 log depth to account for the TL up to range = depth
TLsc2(:,(Close+1):end) = TLsc2(:,(Close+1):end)+ repmat(10*logRmat2(:,Close),1,(length(R)-Close)); 

figure;plot(log10(F),TLc(:,end),log10(F),TLsc1(:,end),log10(F),TLsc2(:,end),log10(F),TLs(:,end),'linewidth',2);
legend('Cylindrical Spreading','Deep case','Shallow case','Spherical Spreading'); 
title('40 kyds');xlabel('Frequency');ylabel('TL');
figure;plot(log10(F),TLc(:,1000),log10(F),TLsc1(:,1000),log10(F),TLsc2(:,1000),log10(F),TLs(:,1000),'linewidth',2);
legend('Cylindrical Spreading','Deep case','Shallow case','Spherical Spreading'); 
title('10 kyds');xlabel('Frequency');ylabel('TL');

