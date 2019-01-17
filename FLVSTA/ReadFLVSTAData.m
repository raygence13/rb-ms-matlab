function  ReadFLVSTAData

%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
headersize = 104;
Fs         = 12500;
Nchans     = 384;
xsecs      = 1;
nPts       = (xsecs*Fs);  % Set up to read x secs of data

time  = ((1:nPts)*(1/Fs));
tfft  = 2^nextpow2(2*Fs);
f     = Fs*(0:tfft/2-1)/tfft;
bins  = 1:tfft/2;

idx_z = 1:4:Nchans/2; % 48 fwd VS elements
idx_y = idx_z + 1;
idx_x = idx_y + 1;
idx_h = idx_x + 1;
%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% File Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fname, pname] = uigetfile('*.rst; ', 'select input RST data file');
if(~ischar(fname))
  return;
end
[nullpath filename theext]=fileparts(fname);
cd(pname);
input_fullname = [pname fname];

fidin    = fopen(input_fullname,'r','l');

fseek(fidin,0,'eof');
filebytes       = ftell(fidin);
timeSampsAvail  = floor((filebytes - headersize) / (4 * Nchans) ); 
fseek(fidin,0,'bof');
fseek(fidin,headersize,'cof');


numBlocks = floor(timeSampsAvail/nPts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% File Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure(1);
set(h,'position',[1 41 1920 947]);

for block = 1:numBlocks
    
    data = fread(fidin,[Nchans,nPts],'float');
    data = data - repmat(mean(data,2),1,nPts); % Remove Mean
    
    fftdata = fft(data,tfft,2)./tfft;
    fftdata = fftdata(:,bins);
    PSDData = fftdata.*conj(fftdata).*2;
    
    h  = data(idx_h,:); % H channels
    x  = data(idx_x,:); % X channels
    y  = data(idx_y,:); % Y channels
    z  = data(idx_z,:); % Z channels
    
    subplot(421)
    imagesc(time,1:length(idx_h),h)
    title('H')
    ylabel('Channel')
    xlabel('time (Secs)');
    
    subplot(422)
    semilogx(f,10*log10(PSDData(idx_h,:)'))
    grid on
    ylim([-160 0]);
    xlim([0 Fs/2]);
    title('H')
    ylabel('Channel')
    xlabel('Frequency (Hz)');
    
    subplot(423)
    imagesc(time,1:length(idx_x),x)
    title('X')
    ylabel('Channel')
    xlabel('time (Secs)');
    
    subplot(424)
    semilogx(f,10*log10(PSDData(idx_x,:)'))
    grid on
    ylim([-160 0]);
    xlim([0 Fs/2]);
    title('X')
    ylabel('Channel')
    xlabel('Frequency (Hz)');
    
    subplot(425)
    imagesc(time,1:length(idx_y),y)
    title('Y')
    ylabel('Channel')
    xlabel('time (Secs)');
    
    subplot(426)
    semilogx(f,10*log10(PSDData(idx_y,:)'))
    grid on
    ylim([-160 0]);
    xlim([0 Fs/2]);
    title('Y')
    ylabel('Channel')
    xlabel('Frequency (Hz)');
    
    subplot(427)
    imagesc(time,1:length(idx_z),z)
    title('Z')
    ylabel('Channel')
    xlabel('time (Secs)');
    
    subplot(428)
    semilogx(f,10*log10(PSDData(idx_z,:)'))
    grid on
    ylim([-160 0]);
    xlim([0 Fs/2]);
    title('Z')
    ylabel('Channel')
    xlabel('Frequency (Hz)');

    pause(0.5)
           
end

fclose(fidin)
return