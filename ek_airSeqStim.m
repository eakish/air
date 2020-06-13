% look at the effect of stimulation on respiratory pressure
% ek 12.21.19

% Load in data
bird = 'or84or24';
date = '121819';
% load(strcat(bird, '_airsac_calls', date, '.mat'))
% load('or84or24_lRA90uA_airsac_callsAndSong121819.mat') % is this still the wrong data?
% load('or84or24_lRA90uA_airsac_stim_song13-Feb-2020.mat')
load('or84or24_lRA90uA_airsac_stim_song121819.mat')
fs = data{1, 1}.syls.Fs;

%% find transitions with stimulation occurring w/in window
% cut out respiratory trace windowing transition - aligned to trans syl offset
trans = 'bb';
win = 500;

for i = 1 : length(trans) - 1 % initialize arrays
    breathArray{i} = [];
    stimArray{i} = [];
end

fs = data{1}.syls.Fs;
% grab 100 ms after syllable offset
airWin = win * fs / 1000;
% x = -win : win * 2 / length(breathArray{j}(:, k)) : win;

for i = 1 : length(data)
%     figure
%     subplot(2, 1, 1)
%     plot(data{i}.song); hold on;
%     subplot(2, 1, 2)
%     plot(data{i}.air);
    if ~isfield(data{i}, 'stimTimes')
        continue
    end
    
    stims = data{i}.stimTimes * 1000 / fs; % in samples
    for j = 1 : length(trans) - 1
    % find transitions
        transInds = strfind(data{i}.syls.labels, strcat(trans(1), trans(j)));
        transTimes = data{i}.syls.offsets(transInds);
        
        % trans times are in samples
        % stim times are in ms
        
        for k = 1 : length(transInds) % getting the wrong indices???
            % check to see if there is a stim occurring in window around transition
            stimTime = intersect(find(stims >= transTimes(k) - win), find(stims <= transTimes(k) + win));
            stimTime = stims(stimTime);
            
            if isempty(stimTime)
                continue
            end
            
%             if length(stimTime) ~= 1
%                 continue
%             end
            
            if transTimes(k) <= win
                continue
            end

            breaths = data{i}.air(transTimes(k) * fs / 1000 - airWin : transTimes(k) * fs / 1000 + airWin);
            breathArray{j} = [breathArray{j} breaths];
            
            stim = (stimTime(1) - transTimes(k));
            stimArray{j} = [stimArray{j} stim];
            
%             figure
%             x = linspace(-win, win, length(breaths));
%             y = breaths;
%             plot(x, y); hold on
%             for s = 1 : length(stimTime)
%                 stim = (stimTime(s) - transTimes(k))
%                 scatter(stim, 0, 100, 'p', 'filled')
%             end
        end           
    end
end
%
transInd = 1;
S = breathArray{transInd}; % need to order by stim time
T = linspace(-win,win,length(S(:, 1)));
stimMat = (stimArray{transInd} + win) * fs / 1000;
% sort by stim onset
[stimMat_sorted, sortInd] = sort(stimMat);
S_sort = S(:, sortInd);

save(strcat(bird, '_', date, '_airsac_sortedStimTime_', trans, '.mat'), 'S_sort', 'stimMat_sorted', 'trans')

%% plot stuff

% figure
% imagesc(S'); hold on
% yStim = (1 : 1 : length(stimArray{transInd}))';
% scatter(stimMat, yStim, 50, 'w', 'filled') 
% title('not sorted')
% 
% xticks = 1 : (length(T) - 1) / 2 : length(T);
% xticklabels = linspace(-win, win, numel(xticks));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% colormap('jet')
% xlabel('Time (ms)')
% ylabel('Trials')


figure
imagesc(S_sort'); hold on
yStim = (1 : 1 : length(stimArray{transInd}))';
scatter(stimMat_sorted, yStim, 50, 'w', 'filled') 
title(['sorted-', trans])

% figure; plot(S)
xticks = 1 : (length(T) - 1) / 2 : length(T);
xticklabels = linspace(-win, win, numel(xticks));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar
xlabel('Time (ms)')
ylabel('Trials Sorted by Stim Time')


%% does respiratory trace look like it was going to another syllable before stim happened?



%% stim triggered average -- what does stimulation do to resp trace on avg?

% broken up by inspiration vs expiration?

% breathing heatmap aligned to insp onset and sorted by stim time


%% respiratory rate vs song tempo
% figure out basal respiratory rate by looking at the number of inspirations over time
% outside of song


% song tempo by looking at syllable rate
sylRate = [];
for i = 1 : length(data)
    % make sure there's only one song per thing
    notSongInds = [strfind(data{i}.syls.labels, 'k')];% strfind(data{i}.syls.labels, '-')];
    notSongInds = sort(notSongInds);
    for j = 2 : length(notSongInds)
        if notSongInds(j - 1) == notSongInds(j)
            continue
        end
        nSyls = length(data{i}.syls.onsets(notSongInds(j - 1) + 1 : notSongInds(j) - 1));
        songDur = (data{i}.syls.onsets(notSongInds(j) - 1) - data{i}.syls.onsets(notSongInds(j - 1) + 1)) / fs; 
        sylRate = [sylRate (nSyls) / songDur]; % syllables per second
    end
end
%%
figure
nhist(sylRate, 'binfactor',20)
ylabel('number songs')
xlabel('syllable Rate')









