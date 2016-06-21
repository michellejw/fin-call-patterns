
clear variables


[pickfiles pickfolder] = uigetfile('*.mseed','Select files','Multiselect','on');

%% If several files are loaded, they're stored as a cell.
% Otherwise it's a character array
if iscell(pickfiles)==1
  lenpickfiles = length(pickfiles)
else
  lenpickfiles = 1;
  pickfiles = {pickfiles};
end

%%

for fdex = 1:length(pickfiles)
  
  [X I] = rdmseed([pickfolder pickfiles{fdex}]);
  
  Net = X(1).NetworkCode;
  Sta = X(1).StationIdentifierCode;
  Chan = X(1).ChannelIdentifier;
  
  datasegmentdex = [1 I.GapBlockIndex' length(X)+1];
  
%   figure(1),clf
%   grid on, hold all
  
  for gdex = 1:length(datasegmentdex)-1
    
    T = cat(1,X(datasegmentdex(gdex):datasegmentdex(gdex+1)-1).t);
    D0 = cat(1,X(datasegmentdex(gdex):datasegmentdex(gdex+1)-1).d);
    D = D0 - mean(D0); 
    
    fs = X(datasegmentdex(gdex)).SampleRate;
    
    % Write wav file. Note: nbits set to 32 (default is 16). rdmseed uses
    % Steim-1 encoding, which I don't entirely understand but I think the 
    % rdmseed.m script says 32-bit encoding (I could be wrong). The reason
    % I even looked into this was because I was getting clipping errors
    % back from wavwrite when using the default 16 bit encoding. After
    % making this change, I stopped getting the clipping errors.
%     wavwrite(D,fs,32,[pickfolder Net '.' Sta '.' Chan '.' datestr(T(1),30) '.wav'])
    timetemp = datevec(T(1)); dsec = num2str(timetemp(6)-floor(timetemp(6)),'%0.3f');
    audiowrite([pickfolder Net '.' Sta '.' Chan '.' datestr(T(1),30) '.' dsec(3:end) '.wav'],D,fs,'BitsPerSample',32);
    
    %plot(T,D)
    
  end
  
%   title(datestr(X(1).t(1),1))
%   datetick 
%   axis tight
%   drawnow
  
  
end