function p= ned(hdg, pit, roll)
% p= ned(hdg, pit, roll)
%

p= [];
close all;

if nargin<1,
 hdg= 0;
 end;
if nargin<2,
 pit= 0;
 end;
if nargin<3,
 roll= 0;
 end;

[hdg, pit, roll]
chdg= cos(hdg*pi/180);
shdg= sin(hdg*pi/180);
cpit= cos(pit*pi/180);
spit= sin(pit*pi/180);
croll= cos(roll*pi/180);
sroll= sin(roll*pi/180);

% rotation vector takes relative bearing vector to vehicle with
% heading, pitch, and roll
% and rotates to true bearing vector
R= [[chdg*cpit, shdg*cpit, -spit]; ...
	[chdg*spit*sroll-shdg*croll, chdg*croll+shdg*spit*sroll, cpit*sroll]; ...
	[chdg*spit*croll+shdg*sroll, -chdg*sroll+shdg*spit*croll, cpit*croll]];
R

az= de= xt= yt= zt= xr= yr= zr= [];

for az0= 0:30:360,
 for de0= 0:0,
  [az0, de0]
  % compute relative bearing
  xt0= cos(az0*pi/180)*cos(de0*pi/180);
  yt0= sin(az0*pi/180)*cos(de0*pi/180);
  zt0= sin(de0*pi/180);
  % rotate to true/earth surface centric
  xyzr= R*[xt0; yt0; zt0];

  az= [az, az0];
  de= [de, de0];
  xt= [xt, xt0];
  yt= [yt, yt0];
  zt= [zt, zt0];
  xr= [xr, xyzr(1)];
  yr= [yr, xyzr(2)];
  zr= [zr, xyzr(3)];
  end;
 end;

figure(1);
plot(xt, yt, 'r+', 'linewidth', 2);

figure(2);
plot3(xt, yt, zt, 'g+', 'linewidth', 2, xr, yr, zr, 'r+', 'linewidth', 2);

hax= gca();
set(hax, 'ydir', 'reverse');
set(hax, 'zdir', 'reverse');

end;

