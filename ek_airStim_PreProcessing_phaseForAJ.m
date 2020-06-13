%% Extracts syllable labels, their on/off times relative to song start, and breathing data
% EK 05/14/19
% EK 12/6/19 - modified to include stim (channels 2 and 3, air - channel 4)
clear all

batchFile = 'batch.txt.dcrd.keep';
% bird = 'rd8';
% date = '110619';
% saveDir = '/home/eszter/Data/airsac/males/rd8/preProcessed/';

bird = 'or84or24';
date = '120619';
experiment = 'lRA50uA';
saveDir = '/home/eszter/Data/airsac/males/or82or24/preProcessed/';
%
% cd '/home/eszter/Data/airsac/males/or82or24/experiment/lRA90uA/2200uF_1ohm_121919/'
cd '/home/eszter/Data/airsac/males/or82or24/experiment/lRA50uA/'
fid=fopen(batchFile,'r');

i = 1;
j = 0;

phase.air = [];
phase.time = [];
phase.callOnsets = [];
phase.callOffsets = [];

figure; hold on
time1 = 0;
while (1)
% for k = 1 : 50
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if ~exist(strcat(fn, '.not.mat'))
        continue;
    end
    j = j + 1;
    [dat,fs]=ReadCbinFile(fn);
    song = dat(:, 1);
    stim1 = dat(:, 2);
    stim2 = dat(:, 3);
    air = dat(:, 4);

    
    syls = open(strcat(fn, '.not.mat')); % figure out first time bin of song file & subtract from on/off times
    if ~isempty(strfind(syls.labels, 'c')) %|| ~isempty(strfind(syls.labels, 'b'))
       
        % use PJ method to find amplitude envelopes
        SPTH = 2000; % toggle to change spectral thresholding
        N = 1024; OVERLAP = 1020; F_low = 500; F_high = 10000; w=hamming(N);
        filtSong = pj_bandpass(song,fs,500,10000,'butterworth');
        rmsSong = rms(filtSong, 2);
%         hold on; plot(filtSong)
        win =  2.5 * fs / 1000;
        envSong = smooth(rmsSong', win);
        
        % filter air sac pressure
%         breathFilt = highpass(air, 1, fs);
%         breathFilt = lowpass(breathFilt, 200, fs);
x = (1:length(air))./fs;
breathFilt = lowpass(air, 100, fs);
% breathFilt1 = highpass(breathFilt, 1, fs);
% breathFilt = highpass(breathFilt, 0.5, fs);
% breathFilt02 = highpass(breathFilt, 0.2, fs);

                

        % overlay traces to look for delay
%         figure
%         plot(x,envSong* 50 - 500)
%         hold on
%         plot(x,breathFilt)


%         x = 1 : 1 : length(air);
%         y = zeros(length(air), 1);
%         plot(x, y, 'k')
%         title(fn)
% %         legend('song', 'air no filter', 'air filtered 1' ,'filterred 0.5','filterred 0.2')
%             legend('song', 'no filter','low pass filtered at 100Hz')

    %   == for stim ===
        index = find(stim2>2000);
        if isempty(index)
            continue
        end
        
        d = diff(index);
        newInd = index(1);

        for l = 1:length(index)-1
            if d(l)>1
                newInd = [newInd index(l+1)];
            end
        end
        
        data{j}.stimTimes = newInd;
        
        data{j}.song = song;
        data{j}.air = breathFilt;
        data{j}.syls = syls;
        
        
        phase.air = [phase.air; breathFilt];
        
        % === find absolute time from file...===
        % datenum(Y,M,D,H,MN,S) 
        date = strsplit(fn, '_');

        day = date{2};
        year = 2000 + str2num(day(end - 1 : end));
        month = str2num(day(3 : 4));
        day = str2num(day(1 : 2));
        
        time = erase(date{3}, '.cbin');
        hour = str2num(time(1 : 2));
        minutes = str2num(time(3 : 4));
        seconds = str2num(time(5 : end));
        
        fileTime = datenum([year, month, day, hour, minutes, seconds]);
        if time1 == 0
            time1 = 1;
%             startTime = fileTime;
            fileStartTime = 0;
        else
            fileStartTime = phase.time(end);
        end
%         fileTime = fileTime - startTime; % calculate time of other files relative to first file time
%         
%         scatter(j, fileTime)
        
        % time relative to file onset
        t = (1000 / fs : 1000 / fs : length(breathFilt) * 1000 / fs) + fileStartTime;
%         pTime = datenum([year .* ones(length(t), 1), month .* ones(length(t), 1), day .* ones(length(t), 1), hour .* ones(length(t), 1), minutes .* ones(length(t), 1), seconds + t']) - startTime;
        phase.time = [phase.time t]; %(t ./ 1000 + fileTime)]; % concatenate
        
        % call onsets - from 'c' -- need to be relative to file time...
        callOnInd = union(strfind(syls.labels, 'ck'), strfind(syls.labels, 'cj'));
        callOnTime = syls.onsets(callOnInd); % add absolute time to here
        phase.callOnsets = [phase.callOnsets; (callOnTime) + fileStartTime]; %; ./ 1000 + fileTime)'];
        
        % call offsets - from 'k' or 'j'
        callOffInd = union(strfind(syls.labels, 'k'), strfind(syls.labels, 'j'));
        callOffTime = syls.offsets(callOffInd); % add absolute time to here...
        phase.callOffsets = [phase.callOffsets; (callOffTime) + fileStartTime]; %; ./ 1000 + fileTime)'];
        
        
        % === sanity check the timing of things ===

        
%         % === plot song, breathing, and syllables and stim!===
%         [S, F, T] = spectrogram(filtSong,w,OVERLAP,N,fs); % for spectrogram
%         T = linspace(0,max(T) * 1000,length(T)) .* fs / 1000;
%         pp = find(abs(S)<=SPTH); %Find entries with very low power + scale up. Makes display nicer.
%         
%         S(pp) = SPTH; S = log(abs(S));
%        
%         figure
%         ax(1) = subplot(2, 1, 1);
%     %         spectrogram(abs(song).^2)
% %         plot(song)
%         imagesc(T, F, S)
%         set(gca,'YDir','normal')
%         ylim([0 10000])
%         colormap(jet);
%         caxis([8 15])
%         title(fn)
%         ax(2) = subplot(2, 1, 2);
%         plot(breathFilt, 'k') % add in a zero line
%         hold on
%         % == to plot stimulations ==
%         scatter(newInd, zeros(length(newInd), 1), 100, 'p', 'r', 'filled')
%         title('airsac')
%         ticks = syls.onsets .* fs / 1000;
%         xticks(ticks) 
%         xLabels = {};
%         for k = 1 : length(syls.labels)
%             xLabels{k} = syls.labels(k);
%         end
%         xticklabels(xLabels)
%         set(gca, 'tickdir', 'out', 'fontsize', 12)
%         linkaxes(ax, 'x')
    end
end
phase.time = phase.time';
figure; plot(phase.time)

fclose(fid)
%
cd(saveDir);
% save(strcat(bird, '_', experiment, '_airsac_callsAndSong', date, '.mat'),'data'); % save entire branchpoint
% save(strcat(bird, '_', experiment, '_airsac_calls', date, '.mat'),'data'); % save entire branchpoint

save(strcat('or84or24_callsForAJ.mat'),'phase'); % save entire branchpoint

















































































