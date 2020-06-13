%% Extracts syllable labels, their on/off times relative to song start, and breathing data
% EK 05/14/19
% EK 12/6/19 - modified to include stim (channels 2 and 3, air - channel 4)
clear all

% === males ===
% batchFile = 'batch.txt';
% experiment = 'airsac';
% bird = 'rd8';
% date = '110719';
% saveDir = '/home/eszter/Data/airsac/males/rd8/preProcessed/';
% cd  '/home/eszter/Data/airsac/males/rd8/calls/'

% batchFile = 'batch.txt.dcrd';
% experiment = 'airsac';
% bird = 'gr98bu75';
% date = '081419';
% saveDir = '/home/eszter/Data/airsac/males/gr98bu75/preProcessed/';
% cd  '/home/eszter/Data/airsac/males/gr98bu75/capacitor/'

bird = 'or84or24';
% 
batchFile = 'batch.txt.dcrd';
date = '120619';
experiment = 'lRA50uA';
cd '/home/eszter/Data/airsac/males/or82or24/experiment/lRA50uA/' % for calls
saveDir = '/home/eszter/Data/airsac/males/or82or24/preProcessed/';

% 
% batchFile = 'batch.txt.dcrd';
% date = '120419';
% experiment = 'baseline';
% cd '/home/eszter/Data/airsac/males/or82or24/baselineBreathing/'  % for baseline and song
% 
% saveDir = '/home/eszter/Data/airsac/males/or82or24/preProcessed/';
% cd '/home/eszter/Data/airsac/males/or82or24/experiment/lRA90uA/2200uF_1ohm_121919/' % for song

% bird = 'rd8';
% batchFile = 'batch.txt';
% date = '110619';
% experiment = 'airsac';
% saveDir = '/home/eszter/Data/airsac/males/rd8/preProcessed/';
% cd '/home/eszter/Data/airsac/males/rd8/calls/'


% % === females ===
% bird = 'or45';
% date = '050819';
% batchFile = 'batch.txt';
% experiment = 'airsac';
% saveDir = '/home/eszter/Data/airsac/FEMALES/or45/preProcessed/';
% cd '/home/eszter/Data/airsac/FEMALES/or45/080319/'


fid=fopen(batchFile,'r');

i = 1;
j = 0;
fs = 32000;
while (1)
% for k = 1 : 50
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    [dat,fs]=ReadCbinFile(fn);
    if isempty(dat)
        continue
    end
    song = dat(:, 1);
    air = dat(:, 2); % for rd8 and before
    
    
%     figure
%     ax(1) = subplot(2, 1, 1);
% %     plot(data{i}.song)
%     N = 1024; OVERLAP = 1020; sigma = 3; F_low = 500; F_high = 10000; filter_type = 'hanningffir';
%     figs = 1;
%     x = 4;
%     y=6;
%     w=hamming(N);
%     SPTH = 300;
%     filtsong=pj_bandpass(song,fs,500,10000,'butterworth');%'hanningffir'
% 
%     [S,F,T] = spectrogram(filtsong,w,OVERLAP,N,fs);
%     T = [1 : length(song)] ./ (fs / 1000);
%     pp = find(abs(S)<=SPTH); %Find entries with very low power + scale up. Makes display nicer.
%     S(pp) = SPTH; S = log(abs(S));
%     hold on
%     imagesc(T,F,S)
%     set(gca, 'ylim', [500 10000], 'xlim', [0, max(T)])
%     colormap jet
%     title(fn)
    
    
%      if ~exist(strcat(fn, '.not.mat')) % ==== get rid of this so it can just run through files
%         continue
%      end
%     
     j = j + 1;
    
%     syls = open(strcat(fn, '.not.mat')); % figure out first time bin of song file & subtract from on/off times
%     if ~isempty(strfind(syls.labels, 'c')) || ~isempty(strfind(syls.labels, 'k')) || ~isempty(strfind(syls.labels, 'i'))
       
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

                

% %         overlay traces to look for delay
%         figure
%         plot(x,envSong* 50 - 500)
%         hold on
%         plot(x,breathFilt)


%         x = 1 : 1 : length(air);
%         y = zeros(length(air), 1);
%         plot(x, y, 'k')
%         title(fn)
%         legend('song', 'air no filter', 'air filtered 1' ,'filterred 0.5','filterred 0.2')
%             legend('song', 'no filter','low pass filtered at 100Hz')
        
        data{j}.song = song;
        data{j}.air = breathFilt;
%         data{j}.syls = syls;

        
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

%     end
end


fclose(fid)
%%
cd(saveDir);
%% save(strcat(bird, '_', experiment, '_airsac_callsAndSong', date, '.mat'),'data'); % save entire branchpoint
save(strcat(bird, '_', experiment, '_airsac_calls', date, '.mat'),'data', '-v7.3'); % save entire branchpoint
% save(strcat(bird, '_', experiment, 'song', date, '.mat'),'data');

















































































