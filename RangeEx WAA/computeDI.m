function DI = computeDI(BP,DE,AZ)
% S. Deeb
% ASSUMES DE/AZ in deg; W/ EQUAL INCREMENTS, & AZ COVERING 360deg
% DE90deg. BP LINEAR abs()^2, NORMALIZED to MAX=1; ITS 1st DIMENSION IS DE

% rd=pi/180;
% de=DE(:)'*rd; az=AZ(:)*rd; 

de=DE(:)'; az=AZ(:); Md = length(DE); Ma=length(AZ);
cd([1,Md]) = 1/3; cd(2:2:Md-1)=4/3; cd(3:2:Md-2)=2/3;
ca([1,Ma]) = 1/3; ca(2:2:Ma-1)=4/3; ca(3:2:Ma-2)=2/3;

ca =(az(2)-az(1))*ca';
cd =(de(2)-de(1))*cd.*cos(de);
DI=db(4*pi/(cd*BP*ca))/2;
