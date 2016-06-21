% run_detectfin_STA.m
% A function to set up parameters for running detectfin.m

clear variables

%% File information

fileprefix = 'NV'; % used for naming the output files
detfolderspecific = 'KEMF_2012-2013_v01/'; % Name for specific results folder
detfolder = ['/burr2/michw/DET_FILES/']; % detections base folder for all results

%% Set up Variables

p.fs = 100;
p.run_id = '20160620a';
p.sta = 'KEMF';

p.resp = 1; % Instrument response
p.tchunk = 1200;
p.tchunkoverlap = 30;
p.wflim = [14 35]; % changed for BB
p.qflim = [7 14];
p.ampfacTable = [0 0.3; 1e11 0.3];
p.mingap = 1;
p.fracindex = 0.1; % changed for BB
p.WINDOW = 125;
p.NOVERLAP = round(0.98*p.WINDOW);
p.NFFT = 125;
p.filttype = 'bandpass';
p.filtorder = 3;
p.SNRthresh = 3;

% Stacking information
p.pickwindow = 200; % window for stacking - seconds
p.dt_up = 13; % time after pick - samples
p.dt_down = 10; % time before pick - samples
p.maxdelay = 50; % max delay time - samples - should be in seconds



%% Make a directory for the detection results files
if exist(detfolder) ~= 7 % if this results directory does not exist yet
    mkdir(detfolder) % make the directory
end

%% Make a directory for the images that are generated during detection
imagefolder = [ detfolder 'images/'];
if exist(imagefolder) ~= 7 % if this images directory does not exist yet
    mkdir(imagefolder) % make the directory
end
p.imagefolder = imagefolder;
p.fileprefix = fileprefix;

%% Load data files: These should be in a specific format: time and
%% amplitude vectors. Time vector should be in matlab serial date format.

[pickfiles pickfolder] = uigetfile('/*.wav','Select hydrophone/OBH/OBS files','Multiselect','on');

% Make sure files are in a cell array
if iscell(pickfiles)==1
    lenpickfiles = length(pickfiles)
else
    lenpickfiles = 1;
    pickfiles = {pickfiles};
end

%% Load call templates
% can generate linear chirps for fin calls here:
% /Users/michw/Research/CODE/DETECTION_CODE_V2

[callfiles, callfolder] = uigetfile('../*.mat','Select call template files','Multiselect','on');

% If several files are loaded, they're stored as a cell.
% Otherwise it's a character array
if iscell(callfiles)==1
    lencallfiles = length(callfiles)
else
    lencallfiles = 1;
    callfiles = {callfiles};
end

for cdex = 1:length(callfiles)
    
    calltemplate(cdex) = load([callfolder callfiles{cdex}]);
    
end

p.calltemplate = calltemplate;

%% Create output text file.
outfilename = [detfolder detfolderspecific(1:end-1) '.run' p.run_id '.txt'];
if exist(outfilename) ~= 2 % if this output file does not exist yet, then create it
    fileID = fopen(outfilename,'w'); % open the file
    fprintf(fileID,'%33s, %10s, %10s, %10s, %10s, %10s\n', ...
        'dettime',...
        'frequency','snr','siglevel','station','run_id');
else
    display('Warning: File already exists. Are you sure you want to append?')
    pause
    fileID = fopen(outfilename,'a'); % open the file for writing, append to existing contents
end

%%
for fdex = 1:length(pickfiles)
    %for fdex = 110:length(pickfiles)
    [data, p.fs] = audioread([pickfolder pickfiles{fdex}]);
    sdex = regexp(pickfiles{fdex},'T')-8; % index at start of date vector
    tfirstsamp = datenum([ str2num(pickfiles{fdex}(sdex:sdex + 3)) ... % Year
        str2num(pickfiles{fdex}(sdex + 4:sdex + 5)) ... % Month
        str2num(pickfiles{fdex}(sdex + 6:sdex + 7)) ... % Day
        str2num(pickfiles{fdex}(sdex + 9:sdex + 10)) ... % Hour
        str2num(pickfiles{fdex}(sdex + 11:sdex + 12)) ... % Min
        (str2num(pickfiles{fdex}(sdex + 13:sdex + 14)) ... 
         + str2num(pickfiles{fdex}(sdex + 15:end-4))) ... % Decimal seconds
        ]);
    
    [dettime, snr, siglev, fmean] = detectfin(data,tfirstsamp,p,0);

    
    for rdex = 1:length(dettime)
        
        dettime_vec = datevec(dettime(rdex));
        yy = dettime_vec(1); mm = dettime_vec(2); dd = dettime_vec(3);
        HH = dettime_vec(4); MM = dettime_vec(5); SS = dettime_vec(6);
        
        % Convert detection time to timestamp format
        dettime_str = [num2str(yy) '-' num2str(mm,'%0.2i') '-' num2str(dd,'%0.2i') ' '...
            num2str(HH,'%0.2i') ':' num2str(MM,'%0.2i') ':' num2str(SS,'%09.6f')];
        
        fprintf(fileID,'%30s, %4.3f, %4.3f, %4.3f, %5s, %12s \n',...
            dettime_str,fmean(rdex),...
            snr(rdex),siglev(rdex),p.sta,p.run_id);
        
    end
    
    
end


fclose(fileID);


