function [ps]=specgram(ts, sr, fn, tod)
% [ps]=specgram(ts, sr, [fn, [tod]])
%
% calculate spectrogram of time series ts at sample rate sr
% when fn present, plot on figure fn
% when tod present, label xaxis as date/time


if (length(ts)<1000),
 avg= 1;
 bs= 32;
elseif (length(ts)<10000),
 avg= 1;
 bs= 128;
elseif (length(ts)<100000),
 avg= 1;
 bs= 512;
else,
 avg= 2;
 bs= 1024;
 end;
step= bs/2;
h= hanning(bs);

ps= [];

% adjust number of steps averaged so specgram image not too large
avgmultiple= floor((length(ts)/step)/10000)+1;
avg= avgmultiple*avg;

if isreal(ts),

avgcount= 0;
for i=1:step:(length(ts)-(bs-1)),
   t= fft(h.*ts(i:(i+bs-1)));
   tt= t.*conj(t);
   if avgcount==0,
      avgtt= tt;
   else
      avgtt= avgtt+tt;
      end;
   avgcount= avgcount+1;
   if avgcount==avg,
      ps= [ps avgtt(1:(bs/2+1))];
      avgcount= 0;
      end;
   end;

if nargin>2,
    figure(fn);
else,
    figure();
    end;

flabel= ((0:(bs/2))/(bs/sr));
if nargin>3 && size(ps,2)*avg*step/sr>600,
    time0= tod;
    day0= floor((tod/(24*3600))+datenum('1/1/1970'));
    tlabel= (1:size(ps,2)-.5)*avg*step/sr+time0;
    tlabel= (tlabel/(24*3600))+datenum('1/1/1970')-day0;
else,
    tlabel= (1:size(ps,2)-1)*avg*step/sr;
    end;
imagesc(tlabel, flabel, 10*log10(ps));
if nargin>3,
    datetick('x','HH:MM:SS');
    end;
axis xy;

else,

avgcount= 0;
for i=1:step:(length(ts)-(bs-1)),
   t= fft(h.*ts(i:(i+bs-1)));
   tt= t.*conj(t);
   tt= [tt((bs/2+1):bs); tt(1:(bs/2))];
   if avgcount==0,
      avgtt= tt;
   else
      avgtt= avgtt+tt;
      end;
   avgcount= avgcount+1;
   if avgcount==avg,
      ps= [ps avgtt(1:bs)];
      avgcount= 0;
      end;
   end;


if nargin>2,
    figure(fn);
else,
    figure();
    end;

flabel= (((-bs/2):(bs/2-1))/(bs/sr));
if nargin>3,
    time0= tod;
    day0= floor((tod/(24*3600))+datenum('1/1/1970'));
    tlabel= (1:size(ps,2)-.5)*avg*step/sr+time0;
    tlabel= (tlabel/(24*3600))+datenum('1/1/1970')-day0;
else,
    tlabel= ((1:size(ps,2))-1)*avg*step/sr;
    end;
imagesc(tlabel, flabel, 10*log10(ps));
if nargin>3,
    datetick('x','HH:MM:SS');
    end;
axis xy;

end;


A= gca();
LW= findall(A, '-property', 'linewidth');
set(LW, 'linewidth', 2);
FW= findall(A, '-property', 'fontweight');
set(FW, 'fontweight', 'bold');
FS= findall(A, '-property', 'fontsize');
set(FS, 'fontsize', 10);

pause(0);
