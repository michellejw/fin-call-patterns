function [dettime, snr, siglev, fmean] = detectfin_v04(data,tfirstsamp,p, makeplotsflag)
%% Function for detection fin whale calls
% Set up parameters and call this function from run_detectfin_*STA*.m
% where *STA* should be replaced with a name that indicates station, date,
% or some identifying notation.
%
% author: Michelle Weirathmueller
% created: 18 August 2014
% modified: 28 August 2014

%% plotting


if makeplotsflag == 1
    h1 = figure(40); clf(h1)
    set(h1,'Visible','on')
    m = ['o','.','*','x','s','^','v'];
    c = ['k','g','b','c','m','y','k'];
end
    


%% Set up call template variables
calltemplate = p.calltemplate;

% Setting up output parameters.  Pre-allocating memory for improved efficiency.
presize = 10000;

dettime = nan(1,presize);
siglev = nan(1,presize);
snr = nan(1,presize);
amplitude = nan(1,presize);
fmean = nan(1,presize);
station = p.sta;

adex = 1;

fs = p.fs;
times = 0:1/fs:length(data)/fs; % time vector in seconds
calltemplate = p.calltemplate;

imagefolder = p.imagefolder; 
fileprefix = p.fileprefix;

tStart = times(1); 
tEnd = times(end); 
nloops = ceil((tEnd-tStart)/(p.tchunk));
looplen = p.tchunk*fs;
localstart = 1;
localend = looplen;

%% Loop through each segment of length tchunk
for ldex = 1:nloops
    if localend < length(data)
        data0 = data(localstart:localend);
        times0 = times(localstart:localend);
        localstart = localend+1;
        localend = localstart+looplen;
    else
        data0 = data(localstart:length(data));
        times0 = times(localstart:length(data));
    end
    
    if length(data0) < p.WINDOW*10
        continue
    end
    
    data0 = data0/p.resp; % Correct for instrument response
    data1 = data0 - nanmean(data0); % demean the data
    % Deal with NaNs and Infs
    data1(isnan(data1)) = 0;
    data1(isinf(data1)) = 0;
    
    %% Build spectrogram
    [~,F0,T,P0] = spectrogram(data1,hann(p.WINDOW),p.NOVERLAP,p.NFFT,fs);
    T = T+times0(1);
    
    Tseis = 0:1/fs:(length(data1)-1)/fs;
    Tseis = Tseis + times0(1);
    Poriginal = P0;
    
    % Find indices of F-vector (frequency vector) for the whale band:
    idWhale = find(F0 > p.wflim(1) & F0 < p.wflim(2));
    idQuake = find(F0 > p.qflim(1) & F0 < p.qflim(2));
    
    % Keep only the whale band from the spectrogram and frequency vectors
    P0 = P0(idWhale,:);
    F = F0(idWhale);
    
    dBB0 = 20*log10(Poriginal);
    dBB = dBB0;
    dBB_dm = dBB0 - repmat(median(dBB0,2),1,length(T));
    
    % Just whale band
    dBWh = dBB(idWhale,:);
    dBQ = dBB(idQuake,:);
    
    % Mean of P matrix (power) for each time, over frequencies in
    % chosen band:
    PWhale = mean(dBWh);
    PQuake = mean(dBQ);
    
    % Log ratio between pWhale and pQuake:
    shiftpos = min([PQuake PWhale])-10; % set min to a positive value to avoid issues when taking logs
    frac = (log10((PWhale-shiftpos)./(PQuake-shiftpos)));
    
    % define the noise level at a certain percentile
    % eg. plot the data like
    % plot(sort(C_demod1))
    % where is it flat before spiking up at the end?
    pct = 90; % setting the desired percentile value
    Psorted = sort(PWhale);
    psortdex = length(Psorted)*((pct/100)); % index in sorted vector of noise floor
    thresh = Psorted(round(psortdex));
    
    % Set display threshold to 2 standard deviations from mean:
    %thresh = mean(dBWh(:)) + 2*std(dBWh(:));
    % Original method used by Will and Dax:
    dBW = max(dBWh - thresh,0);
    % Create a matrix of frequencies that is the same size as dBW:
    fdB = repmat( F , 1 , size(dBW,2) );
    % Find frequency mean and standard deviation
    fmean0 = sum(fdB.*dBW) ./ sum(dBW);
    fstd = sqrt( sum(dBW.*(fdB-repmat(fmean0, length(idWhale),1)).^2) ./sum(dBW));
    
    % For plotting, get percentiles 
    alldbb = dBB(:);
    pct_dBB = 5;
    dBBsorted = sort(alldbb);
    dBBsortdex = length(dBBsorted)*(pct_dBB/100);
    dBB_pctval = dBBsorted(round(dBBsortdex));
    
    
    %% Filter data (used for plotting only, not for detection)
    % Normalized cutoff frequency (cut-off frequency/(0.5*sample rate).
    % Wn must be less than 1.  If it is not, an error message is generated,
    % and the program exits.
    wn = 2 * p.wflim/fs;
    if any(find(wn > 1))
        error(['The filter cut-off frequency is too high for a' ...
            ' sample rate of ' num2str(fs)]);
    else
        [b,a] = butter( p.filtorder, wn, p.filttype );
        datafilt = filtfilt(b, a, data1);
    end
    
    
    %% DETECT CALLS
        
%        load([callfolder callfiles{tdex}])
        
        fs_template = calltemplate.fs_template;
        tvec = calltemplate.tvec;
        fincall = calltemplate.fincall;
        cf = calltemplate.cf;
        bw = calltemplate.bw;
        chlength = calltemplate.chlength;
        
        % Ensure sample rate is the same for template and data
        if fs ~= fs_template
            error(['Template sample rate is not the same as data sample rate'])
        end
        
        
        %% Matched filter
        
        % Matched filtering routine
        [B,A] = butter( p.filtorder , (bw*1.5)*2/fs);
        tx1 = cos(2*pi*cf*tvec).*fincall;
        tx2 = sin(2*pi*cf*tvec).*fincall;
        I_tx = filtfilt(B,A,tx1);
        Q_tx = filtfilt(B,A,tx2);
        tx_demod = I_tx + i*Q_tx;
        
        s1 = cos(2*pi*cf*Tseis).*data1';
        s2 = sin(2*pi*cf*Tseis).*data1';
        I_sig = filtfilt(B,A,s1);
        Q_sig = filtfilt(B,A,s2);
        sig_demod = I_sig + i*Q_sig;
        
        % cross correlate
        C_demod1 = abs(xcorr(datafilt,fincall))/fs;
        C_demod1 = C_demod1(ceil(length(C_demod1)/2):end);
        C_demod1(1:length(fincall)) = 0;
        try
            C_demod1(end-length(fincall):end) = 0;
        catch
            display('stop')
        end
        
        % Find indices in amp vector where a given element is greater than the
        % surrounding elemenst on either side (very local maxima):
        Cdex = (find( C_demod1(2:end-1) > C_demod1(1:end-2) & C_demod1(2:end-1) > C_demod1(3:end)))+1;
        
        tWhalem0 = Tseis(Cdex);
        aWhalem0 = C_demod1(Cdex);
        
        % do fractional analysis:  compare with earthquake amplitude.
        % interpolate to find frac values at xcorr detection locations:
        % fracdelay = chlength/2;
        % fracdelay = 1;
        fracdelay = p.calltemplate.chlength;
        frac2 = interp1(T-fracdelay,frac,...
            tWhalem0,'linear','extrap');
        isWhale = find(abs(frac2 > p.fracindex));
        tWhalem = tWhalem0(isWhale);
        aWhalem = aWhalem0(isWhale);
        
        
        biggap = 10;
        % test for minimum gaps (pickpeaks.m) - careful - this is also a native
        % Matlab function
        [TBigm detBigm] = pickpeaks_mingap(tWhalem, aWhalem, biggap);
        
        [sortT sortdex] = sort(TBigm,'ascend');
        TBig0 = sortT;
        detBig0 = detBigm(sortdex);
        
        
        clear SNRbig
        
        if length(TBig0) == 0
            SNRbig = [];
            noisechunk = [];
            noisefloor = [];
        else
            
            % define the noise level at a certain percentile
            % eg. plot the data like
            % plot(sort(C_demod1))
            % where is it flat before spiking up at the end?
            pct = 90; % setting the desired percentile value
            Csorted = sort(C_demod1);
            nsdex = length(Csorted)*((pct/100)); % index in sorted vector of noise floor
            noisefloor = Csorted(round(nsdex));
            
            
            % determine SNR for each pick
            
            SNRbig = 20*log10(detBig0/noisefloor); % Will says it should be 20, because output of xcorr is not proportional to power.
            
            
            
        end
        
%         SNRgood = find(SNRbig > p.SNRthresh);
%         TBig = TBig0(SNRgood);
%         detBig = detBig0(SNRgood);
        
        %% Extract RMS signal and noise & calculate improved frequency
        % estimate
        sigdB = nan(size(TBig0));
        noisedB = nan(size(TBig0));
        fmean_final = nan(size(TBig0));
        newTBig0 = nan(size(TBig0));
        
%         figure(3),clf
%         plot(Tseis,datafilt)
%         hold on, grid on
        
        for rdex = 1:length(TBig0)
            
            thiscall = find(Tseis > TBig0(rdex)-(bw/2) & Tseis < TBig0(rdex) + (bw/2));
            prenoise = find(Tseis > TBig0(rdex)-(bw) & Tseis < TBig0(rdex)-(bw/2));
            
            % Find the maximum amplitude within the call window
            [~, maxampsubdex] = max(datafilt(thiscall));
            newTBig0(rdex) = Tseis(thiscall(maxampsubdex));
            % secondary subset (we need to do this because I'm using a
            % longer chirp length to capture different frequency calls.
            nsamps = 0.5 * fs;
            thiscall = thiscall(thiscall > thiscall(maxampsubdex)-nsamps ...
                & thiscall < thiscall(maxampsubdex)+nsamps);
            
            sigrms = sqrt(mean(datafilt(thiscall).^2));
            noiserms = sqrt(mean(datafilt(prenoise).^2));
            
            sigdB(rdex) = 20*log10( sigrms );
            noisedB(rdex) = 20*log10( noiserms );
            
            %
%             figure(3)s
%             plot(Tseis(thiscall),datafilt(thiscall),'r')

            
            % Extract a subset of the spectrogram around the detection,
            % using the extents of the time and frequeny of the matched
            % filter chirp.
%           % Since the detection time corresponds to the onset time, we
%           begin the sub-window at that time and go for the length of the
%           chirp. We extend slightly beyond the actual call bandwidth in
%           case the call is slightly skewed, or has a broader bandwidth
%           than the model chirp.
            thiscallT = T >= newTBig0(rdex) & T <= newTBig0(rdex) + (chlength*2);
            thiscallF = F0 >= cf-(bw/2) & F0 <= cf + (bw/2);
            spectrogram_subset = dBB_dm(thiscallF,thiscallT);
            spectrogram_subset(spectrogram_subset<0) = 0;
            
%             figure(11),clf
%             imagesc(T2(thiscallT),F0(thiscallF),spectrogram_subset)
%             datetick
            
            % Now calculate the weighted mean frequency for each time slice
            % in the subset spectrogram
%             fmat = repmat(F0(thiscallF),1,length(thiscallT));
%             fmean_sub = sum(fmat.*spectrogram_subset) ./ sum(spectrogram_subset);
            fmat = F0(thiscallF);
            weights =  sum(spectrogram_subset,2)-min( sum(spectrogram_subset,2));
            fmean_sub = sum(fmat.*weights) / sum(weights);
            
            fmean_final(rdex) = median(fmean_sub);
            
        end
        
        % Re do the min-gap search since timing may have shifted in
        % previous step
        [newTBig0_mingap, ~] = pickpeaks_mingap(newTBig0, SNRbig, biggap);
        % Find indices for detections that remain.
        [~,mingapdex,~] = intersect(newTBig0,newTBig0_mingap); 
        SNRbig = sigdB - noisedB;
        SNRgooddex = find(SNRbig > p.SNRthresh);
        SNRgood = SNRbig(intersect(SNRgooddex,mingapdex));
        TBig = newTBig0(intersect(SNRgooddex,mingapdex));
        detBig = detBig0(intersect(SNRgooddex,mingapdex));

        sigdB = sigdB(intersect(SNRgooddex,mingapdex));
        noisedB = noisedB(intersect(SNRgooddex,mingapdex));
        fmean_final = fmean_final(intersect(SNRgooddex,mingapdex));
        
        
        %%
        if adex + (length(TBig-1)) <= length(dettime)
            
            fmean(adex : adex + (length(TBig)-1)) = fmean_final;
            siglev(adex : adex + (length(TBig)-1)) = sigdB;
            dettime(adex : adex + (length(TBig)-1)) = (TBig/86400) + tfirstsamp;
            snr(adex : adex + (length(TBig)-1)) = SNRgood;
            

            adex = adex + length(TBig);
            
        else
            
            fmean = [fmean nan(1,presize)];
            siglev = [siglev nan(1,presize)];
            dettime = [dettime nan(1,presize)];
            snr = [snr nan(1,presize)];
            
            fmean(adex : adex + (length(TBig)-1)) = fmean_final;
            siglev(adex : adex + (length(TBig)-1)) = sigdB;
            dettime(adex : adex + (length(TBig)-1)) = (TBig/86400) + tfirstsamp;
            snr(adex : adex + (length(TBig)-1)) = SNRgood;

            adex = adex + length(TBig);
            
        end
        
        
        if makeplotsflag == 1
            %%
            clf(h1)
            ax(1) = subplot(411)
            imagesc(T/86400 + tfirstsamp,F0,dBB),axis xy
            hold on
            plot(TBig/86400 + tfirstsamp,fmean_final,[m(1) c(1)])
            datetick
            axis tight
            if fs <= 50
                ylim([0 25])
            else
                ylim([0 50])
            end
            ylabel('Frequency (Hz)','fontsize',14)
            title(datestr(tfirstsamp+(times0(1)/86400),31))
            cax_auto = caxis;
            caxis([dBB_pctval,cax_auto(2)])
            %     colorbar
            
            ax(2) = subplot(412)
%             plot(Tseis/86400+tfirstsamp,C_demod1), hold on
%             plot(TBig/86400 + tfirstsamp, detBig,'ok')
            plot(Tseis/86400+tfirstsamp,datafilt), hold on
            plot(TBig/86400+tfirstsamp,interp1(Tseis/86400+tfirstsamp,...
                datafilt,TBig/86400+tfirstsamp),'ok')
            grid on
            datetick
            axis tight
            
            ax(3) = subplot(413)
            plot(TBig/86400 + tfirstsamp,SNRgood,[m(1) c(1)])
            datetick
            hold on, grid on,
            
            ax(4) = subplot(414)
            %plot(T/86400 + tfirstsamp-(fracdelay/86400),frac), grid
            plot(T/86400 + tfirstsamp,frac), grid
            datetick
            
            linkaxes(ax,'x')
            
%            saveas(h1,[imagefolder fileprefix datestr(tfirstsamp+(times0(1)/86400),30)],'png')
        end
        %%
        
        clear fmean_final fstd_final SNRbig SNRgood TBig Cdemod1 T2 detBig noisedB
        
    

    
end


display('stop')

nandex = isnan(dettime);
dettime(nandex) = [];
snr(nandex) = [];
siglev(nandex) = [];
fmean(nandex) = [];





% lowsnrdex = snr<5;
% dettime(lowsnrdex) = [];
% snr(lowsnrdex) = [];
% siglev(lowsnrdex) = [];
% fmean(lowsnrdex) = [];

%ipi = [0 diff(dettime*86400)];
