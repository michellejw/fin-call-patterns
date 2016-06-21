% Script for downloading large amounts of data from the Iris website
% (a month at a time from a particular instrument)
%
% For examples, see:
% http://service.iris.edu/irisws/timeseries/docs/1/builder/
%
% for more info on Cascadia stations, see:
% http://www.iris.washington.edu/gmap/_CASCADIA_OBS

clear variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Fill in the following information:          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

network = 'NV' % 7D is CIET network, NV is neptune canada/ONC

stations{1} = 'KEMF'; % Station name (can be a cell 
% array of stations if channel names do not change)

% Desired channels: these may change - check station 
% info on website!
channels{1} = 'EHZ'; 
%channels{2} = 'BH1';
%channels{3} = 'BH2';

% Start and end dates, converted to matlab datenum 
% format
startdate = datenum('01-Oct-2012','dd-mmm-yyyy');
enddate = datenum('01-Apr-2013','dd-mmm-yyyy');

irisformat = 'miniseed'; % miniseed, audio (wav), ascii
file_extension = '.mseed'; % file extension that matches 
% your desiredformat (e.g. '.wav', '.mseed', '.txt', etc)

% irisformat = 'miniseed'; % miniseed, audio (wav), ascii
% file_extension = '.mseed'; % file extension that matches 
% % your desiredformat (e.g. '.wav', '.mseed', '.txt', etc)

% fs = 125; % audio sample rate - used when downloading 
% .wav files - find this on the iris website for specific
% instruments. Not sure why Iris download tool can't 
% populate this field automatically, but it may be that 
% sample rate can change on a given station/instrument, and 
% you can't have a variable sample rate when using a .wav
% file.

% Download size: number of days per file: maximum = 30 
% (max size imposed by Iris).
dayincr = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if enddate < startdate
  error('Start data is after end date')
end

if dayincr > 30
  error('Day increment is too large')
end

localstart = startdate;
localend = localstart+dayincr;

while localend <= enddate % loop through data 30 days at a time
  
  for stadex = 1:length(stations) % loop through stations
    
    for chandex = 1:length(channels) % loop through channels
      
      % build URL for particular station/channel/time period. There are two
      % lines below, one that works for WAV audio files, and one that
      % works for miniseed files. I need to make this automatic.
      
      % WAV files:
%       url = ['http://service.iris.edu/irisws/timeseries/1/query?net=' network '&sta=' stations{stadex} '&loc=--&cha=' channels{chandex} '&start=' datestr(localstart,29) ...
%         'T00:00:01&end=' datestr(localend,29) 'T00:00:01&output=' irisformat '&audiosamplerate=' num2str(fs) ];

%       % miniseed files:
      url = ['http://service.iris.edu/irisws/timeseries/1/query?net=' network '&sta=' stations{stadex} '&loc=--&cha=' channels{chandex} '&start=' datestr(localstart,29) ...
        'T00:00:01&end=' datestr(localend,29) 'T00:00:01&output=' irisformat ];
      
      % download data to specified file
      try
        urlwrite(url, [stations{stadex} channels{chandex} '.' datestr(localstart,30) file_extension]);
      catch
        continue
      end
      
    end
    
  end
  
  localstart = localend;
  localend = localend + dayincr;
  
end

% download remaining fraction of month

localend = enddate;

if localstart ~= localend
  % WAV file:
%   url = ['http://service.iris.edu/irisws/timeseries/1/query?net=' network '&sta=' stations{stadex} '&loc=--&cha=' channels{chandex} '&start=' datestr(localstart,29) ...
%     'T00:00:01&end=' datestr(localend,29) 'T00:00:01&output=' irisformat '&audiosamplerate=' num2str(fs)];
  
  % % Miniseed file:
  url = ['http://service.iris.edu/irisws/timeseries/1/query?net=' network '&sta=' stations{stadex} '&loc=--&cha=' channels{chandex} '&start=' datestr(localstart,29) ...
    'T00:00:01&end=' datestr(localend,29) 'T00:00:01&output=' irisformat ];
  
  urlwrite(url, [stations{stadex} channels{chandex} '.' datestr(localstart,30) file_extension]);
  
  
end


