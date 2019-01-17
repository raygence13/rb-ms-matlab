function readFLVSTARoll
rollidx = 4:4:192;

[fname, pname] = uigetfile('*.dat; ', 'select input RST data file');
if(~ischar(fname))
  return;
end
[nullpath filename theext]=fileparts(fname);
cd(pname);
input_fullname = [pname fname];

fidin    = fopen(input_fullname,'r','l');

Header.HeaderBytes = fread(fidin,1,'int');
Header.DateTime    = fread(fidin,1,'int');
Header.Chans       = fread(fidin,1,'int');
Header.ModId       = fread(fidin,1,'int');

RollBlock  = 13; 

fseek(fidin,0,'eof');
filebytes       = ftell(fidin);
timeSampsAvail  = floor((filebytes - Header.HeaderBytes) / (4 * Header.Chans*4) ); 
fseek(fidin,0,'bof');
fseek(fidin,Header.HeaderBytes,'cof');


numBlocks = floor(timeSampsAvail/RollBlock);

for block = 1:numBlocks
    
    RollData = fread(fidin,[Header.Chans*4 RollBlock],'float');
    
    roll(:,block) = unwrap(mean(RollData(rollidx,:),2)); 
    
end
figure(2)
plot(roll')
grid on

fclose(fidin);

return