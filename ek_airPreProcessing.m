%% Extracts syllable labels, their on/off times relative to song start, and breathing data
% EK 05/14/19

batchFile = 'batch.txt.dcrd.keep';
% bird = 'rd8';
% date = '110619';
% saveDir = '/home/eszter/Data/airsac/males/rd8/preProcessed/';

bird = 'or84or24';
date = '120819';
saveDir = '/home/eszter/Data/airsac/males/or84or24/preProcessed/';
%%
fid=fopen(batchFile,'r');

i = 1;
j = 0;
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    j = j + 1;
    [dat,fs]=ReadCbinFile(fn);
    song = dat(:, 1);
    air = dat(:, 2);
    syls = open(strcat(fn, '.not.mat')); % figure out first time bin of song file & subtract from on/off times
    if ~isempty(strfind(syls.labels, 'k')) %|| ~isempty(strfind(syls.labels, 'b'))
       
        % use PJ method to find amplitude envelopes
        N = 1024; OVERLAP = 1020; sigma = 3; F_low = 500; F_high = 10000; filter_type = 'hanningffir';
        filtSong = highpass(song, F_low, fs);
%         figure; plot(filtSong)
        filtSong = lowpass(filtSong, F_high, fs);
%         hold on; plot(filtSong)
        filtSong = rms(filtSong, 2);
%         hold on; plot(filtSong)
        win =  2.5 * fs / 1000;
        envSong = smooth(filtSong', win);
        
        % filter air sac pressure
        breathFilt = highpass(air, 1, fs);
        breathFilt = lowpass(breathFilt, 200, fs);
        
        % overlay traces to look for delay
        figure
        plot(envSong* 50 - 500)
        hold on
        plot(air)
        hold on
        plot(breathFilt)
        hold on
        x = 1 : 1 : length(song);
        y = zeros(length(song), 1);
        plot(x, y, 'k')
        legend('song', 'air no filter', 'air filtered')
        
        
        data{j}.song = song;
        data{j}.air = breathFilt;
        data{j}.syls = syls;
        
%     % plot song, breathing, and syllables
%     figure
%     ax(1) = subplot(2, 1, 1)
% %         spectrogram(abs(song).^2)
%     plot(song)
%     title(fn)
%     ax(2) = subplot(2, 1, 2)
%     plot(air)
%     title('airsac')
%     ticks = syls.onsets .* fs / 1000;
%     xticks(ticks) 
%     xLabels = {};
%     for k = 1 : length(syls.labels)
%         xLabels{k} = syls.labels(k);
%     end
%     xticklabels(xLabels)
%     set(gca, 'tickdir', 'out', 'fontsize', 12)
%     linkaxes(ax, 'x')
    end
end
    
fclose(fid)
%%
cd(saveDir);
save(strcat(bird, '_airsac_', date, '.mat'),'data'); % save entire branchpoint








