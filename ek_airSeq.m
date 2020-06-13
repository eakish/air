%% get air sac traces for specific transitions

% Load in data
bird = 'rd8';
date = '110619';
load(strcat(bird, '_airsac_', date, '.mat'))
fs = data{1}.syls.Fs;
% load('pk53gr18_051219_calls.mat')

%
% === gr8bu
% trans = 'ajhib';
% trans = 'hij';
% trans = 'bcj'; % stereotyped transitions

% == control for place in song == 
% preTrans = 'hjkabcbjk'; % syllables that should precede transition

% === rd8 ===
trans = 'kc';

% == control for place in song == 
preTrans = ''; % syllables that should precede transition

% == control for gap duration ==
minGap = 0;
maxGap = 1000; 

% plot a single example song and respiratory trace
N = 1024; OVERLAP = 1020; sigma = 3; F_low = 500; F_high = 10000; filter_type = 'hanningffir';
figs = 1;
x = 4;
y=6;
w=hamming(N);
SPTH = 300;
filtsong=pj_bandpass(data{1}.song,fs,500,10000,'butterworth');%'hanningffir'

figure; 
bigFig(1) = subplot(2, 1, 1); 
[S,F,T] = spectrogram(filtsong,w,OVERLAP,N,fs);
T = [1 : length(data{1}.air)] ./ (fs / 1000);
pp = find(abs(S)<=SPTH); %Find entries with very low power + scale up. Makes display nicer.
S(pp) = SPTH; S = log(abs(S));
hold on
imagesc(T,F,S)
set(gca, 'ylim', [500 10000], 'xlim', [-0.5, 0.5])
colormap jet


% plot(data{1}.song); hold on
y = [1 : length(data{1}.air)] ./ (fs / 1000);
bigFig(2) = subplot(2, 1, 2); 
plot(y, data{1}.air, 'k', 'linewidth', 1.5);
%
set(gca, 'tickdir', 'out', 'xlim', [min(y) max(y)])
linkaxes(bigFig, 'x')


%% ===== need to do ==== %%%%%
% control of motif repeat number - done!
% control for gap length - look at gap lengths between 100 & 200 ms? -done!
% use ROC analysis to figure out when respiratory signals diverge
% is there a change in syllable amplitude across song? - plot amplitude at
    % a specific syllable (a?) wrt to syllable index...


%% iterate through all songs and find these transitions
% at these transitions put 200 ms after syllable onset into array
breathArray{1}.songStop.air = []; % for song stoppings
for i = 2 : length(trans) % initialize arrays, where last one is song stopping
%     eval(strcat('breathArray.', trans(i), ' = []'));
    breathArray{i}.sylOffset = [];
    breathArray{i}.sylOnset = [];
    breathArray{i}.gapOffset = [];
    
    breathParams{i}.inspPeakAmp = [];
    breathParams{i}.inspPeakTime = [];
    breathParams{i}.inspVolume = [];
    breathParams{i}.inspDur = [];
    breathParams{i}.gapDur = [];
    
%     sylParams.amplitude = [];
%     sylParams.expiration = [];
%     sylParams.pitch = [];
%     sylParams.sylNum = []; % location in song sequence
end

for i = 1 : 30 %length(data)
    sylParams{i}.amplitude = [];
    sylParams{i}.expiration = [];
    sylParams{i}.pitch = [];
    sylParams{i}.sylNum = []; % location in song sequence
end


% grab windows before and after syllable offset/onset
postWindow = 200 .* (fs / 1000);
preWindow = 100 .* (fs / 1000);

% figure
songStops = [];
% preTrans = 'hjkabcbjk'; % syllables that should precede transition
preTrans = '';
% == control for gap duration ==
minGap = 200;
maxGap = 300; 
for i = 1 : 30 %length(data)
    if isempty(data{i})
        continue
    end
    
%     data{i}.air = data{i}.air + 1000; % for rd8 because weird offset?

%     figure
%     subplot(2, 1, 1)
%     plot(data{i}.song); hold on;
%     subplot(2, 1, 2)
%     plot(data{i}.air);

%     % ====== song stopping!!! =====
%     if strcmp(data{i}.syls.labels(end), trans(1))
%         songStopInd = data{i}.syls.offsets(end) .* (fs / 1000);
%         breathArray{1}.songStop.air = [breathArray{1}.songStop.air data{i}.air(songStopInd - preWindow : songStopInd + postWindow)];
%     end

        % syllable parameters
    sylParams{i}.sylNum = [sylParams{i}.sylNum strfind(data{i}.syls.labels, trans(1))];

    for k = 1 : length(sylParams{i}.sylNum)
        sylN = sylParams{i}.sylNum(k);
        syl = data{i}.song(data{i}.syls.onsets(sylN) * fs/1000 : data{i}.syls.offsets(sylN) * fs/1000);
%         figure; plot(syl)
        exp = data{i}.air(data{i}.syls.onsets(sylN) * fs/1000 : data{i}.syls.offsets(sylN) * fs/1000);
        sylParams{i}.expiration = [sylParams{i}.expiration max(exp)]; 
        sylParams{i}.amplitude = [sylParams{i}.amplitude max(abs(syl))];
%             sylParams{j}.pitch(k) = []; -- use paul's code on egret for this
    end     
%     figure(10); hold on
%     scatter(sylParams{i}.sylNum, sylParams{i}.amplitude)
%     
%     figure(11); hold on
%     scatter(sylParams{i}.amplitude, sylParams{i}.expiration)

    for j = 2 :  length(trans)
    % find transitions
        transInds = strfind(data{i}.syls.labels, strcat([preTrans trans(1)], trans(j))) + length(preTrans);
        transSylOffset = data{i}.syls.offsets(transInds) .* (fs / 1000);
        transSylOnset = data{i}.syls.onsets(transInds) .* (fs / 1000);
        transGapOffset = data{i}.syls.onsets(transInds + 1) .* (fs / 1000);

%         figure
        for k = 1 :  length(transInds) - 1         
            
            gap = data{i}.air(transSylOffset(k) : transGapOffset(k)) .* (fs / 1000);
%             plot(gap); hold on

            % to control for gap duration
            if length(gap) < minGap | gap > maxGap
                continue
            end
            
            insp = gap(find(gap < 0)); % change this
            if isempty(insp)
                continue
            end
            
            breathArray{j}.sylOffset = [breathArray{j}.sylOffset data{i}.air(transSylOffset(k) - preWindow : transSylOffset(k) + postWindow)]; %data{i}.air(transSylOffset(k) - window : transSylOffset(k) + window)]; % also add in song stopping...
            breathArray{j}.sylOnset = [breathArray{j}.sylOnset data{i}.air(transSylOnset(k) - preWindow : transSylOnset(k) + postWindow)];
            breathArray{j}.gapOffset = [breathArray{j}.gapOffset data{i}.air(transGapOffset(k) - preWindow : transGapOffset(k) + postWindow)];
            
            breathParams{j}.gapDur = [breathParams{j}.gapDur length(gap) ./(fs / 1000)]; %data{i}.syls.onsets(transInds(k) + 1) - data{i}.syls.offsets(transInds(k))];
            
            
            
%             x = [1 : length(gap)] ./ (fs / 1000);
%             figure(j - 1)
%             plot(x, gap, 'color', cmap(j - 1, :)); hold on
            
%             plot(insp); hold on;
            
%             [Y, I] = min(gap); %min((insp));
            [Y, I] = findpeaks(-gap,'MinPeakHeight', 0.5e5, 'MinPeakProminence',0.05e5);
            if isempty(Y)
                continue
            end
            Y = -Y(1);
            I = I(1) / (fs / 1000);
            breathParams{j}.inspPeakAmp = [breathParams{j}.inspPeakAmp Y];
            breathParams{j}.inspPeakTime = [breathParams{j}.inspPeakTime I / (fs / 1000)];
           
            breathParams{j}.inspVolume = [breathParams{j}.inspVolume sum(insp)];
            
%             x = (1 : 1 : length(gap)) ./ (fs / 1000);
%             figure
%             plot(x, gap); hold on
%             scatter(I, Y, 100, 'filled')
%             xlabel('Time from gap offset (ms)')
%             ylabel('pressure (mv)')
%             set(gca, 'tickdir', 'out', 'fontsize', 20)
            
            
        end           
    end
end

% figure(10)
% xlabel('Syllable Number in Song')
% ylabel('Syllable Amplitude')
% set(gca, 'tickdir', 'out', 'fontsize', 20)
% 
% figure(11)
% xlabel('Syllable Amplitude')
% ylabel('Expiratory Amplitude')
% set(gca, 'tickdir', 'out', 'fontsize', 20)

%% Make figures!
s = 80; % size of scatter points
cmap = colormap('lines');
for i = 2 : length(trans)
    figure(4) % look at relationship between gap duration and peak inspiratory amplitude
    scatter(breathParams{i}.gapDur, abs(breathParams{i}.inspPeakAmp), s, cmap(i - 1, :), 'linewidth', 2); hold on;
    xlabel('Gap Duration (ms)')
    ylabel('Peak Amplitude of Inspiration')
    
    % remove outliers
    [gapDur_noOutliers, outlierIdx] = rmoutliers(breathParams{i}.gapDur);
    inspPeakAmp_noOutliers = breathParams{i}.inspPeakAmp(find(outlierIdx == 0));
    inspVolume_noOutliers = breathParams{i}.inspVolume(find(outlierIdx == 0));
    inspPeakTime_noOutliers = breathParams{i}.inspPeakTime(find(outlierIdx == 0));
    
    % regression to test for correlation
    stats.linReg.gapDur.inspPeakAmp{i}.r = corr(breathParams{i}.gapDur', breathParams{i}.inspPeakAmp');
    [stats.linReg.gapDur.inspPeakAmp{i}.p, stats.linReg.gapDur.inspPeakAmp{i}.S] = polyfit(gapDur_noOutliers,inspPeakAmp_noOutliers,1);
    [stats.linReg.gapDur.inspPeakAmp{i}.f, stats.linReg.gapDur.inspPeakAmp{i}.delta] = polyval(stats.linReg.gapDur.inspPeakAmp{i}.p, gapDur_noOutliers, stats.linReg.gapDur.inspPeakAmp{i}.S');
    h = plot(gapDur_noOutliers, stats.linReg.gapDur.inspPeakAmp{i}.f, 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspPeakAmp{i}.f + (2*stats.linReg.gapDur.inspPeakAmp{i}.delta), '--', 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspPeakAmp{i}.f - (2*stats.linReg.gapDur.inspPeakAmp{i}.delta), '--', 'color', cmap(i - 1, :)); hold on;
%     label(h, strcat('R^2 = ', num2str(stats.linReg.gapDur.inspPeakAmp{i}.r ^2)), 'location', 'west', 'fontsize', 13, 'fontweight', 'bold')
    
    figure(5) % look at relationship between gap duration and inspiratory volume
    scatter(breathParams{i}.gapDur, abs(breathParams{i}.inspVolume), s, cmap(i - 1, :), 'linewidth', 2); hold on;
    xlabel('Gap Duration (ms)')
    ylabel('Inspiratory Volume') % actually do inspiratory volume instead of determining by syl offset/onset
    
     % regression to test for correlation
    stats.linReg.gapDur.inspVolume{i}.r = corr(breathParams{i}.gapDur', breathParams{i}.inspVolume');
    [stats.linReg.gapDur.inspVolume{i}.p, stats.linReg.gapDur.inspVolume{i}.S] = polyfit(gapDur_noOutliers,inspVolume_noOutliers,1);
    [stats.linReg.gapDur.inspVolume{i}.f, stats.linReg.gapDur.inspVolume{i}.delta]  = polyval(stats.linReg.gapDur.inspVolume{i}.p, gapDur_noOutliers, stats.linReg.gapDur.inspVolume{i}.S');
    h = plot(gapDur_noOutliers, stats.linReg.gapDur.inspVolume{i}.f, 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspVolume{i}.f + 2*stats.linReg.gapDur.inspVolume{i}.delta, '--', 'color', cmap(i - 1, :))
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspVolume{i}.f - 2*stats.linReg.gapDur.inspVolume{i}.delta, '--', 'color', cmap(i - 1, :)); hold on;
%     label(h, strcat('R^2 = ', num2str(stats.linReg.gapDur.inspVolume{i}.r ^2)), 'location', 'west', 'fontsize', 13, 'fontweight', 'bold')
    
    figure(6) % look at relationship between gap duration and peak inspiratory time
    scatter(breathParams{i}.gapDur, breathParams{i}.inspPeakTime, s, cmap(i - 1, :), 'linewidth', 2); hold on;
    xlabel('Gap Duration (ms)')
    ylabel('Inspiratory Peak Time From Syllable Offset (ms)')
    
     % regression to test for correlation
    stats.linReg.gapDur.inspPeakTime{i}.r = corr(breathParams{i}.gapDur', breathParams{i}.inspPeakTime');
    [stats.linReg.gapDur.inspPeakTime{i}.p, stats.linReg.gapDur.inspPeakTime{i}.S] = polyfit(gapDur_noOutliers, inspPeakTime_noOutliers,1);
    [stats.linReg.gapDur.inspPeakTime{i}.f, stats.linReg.gapDur.inspPeakTime{i}.delta] = polyval(stats.linReg.gapDur.inspPeakTime{i}.p, gapDur_noOutliers, stats.linReg.gapDur.inspPeakTime{i}.S');
    h = plot(gapDur_noOutliers, stats.linReg.gapDur.inspPeakTime{i}.f, 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspPeakTime{i}.f + 2*stats.linReg.gapDur.inspPeakTime{i}.delta, '--', 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspPeakTime{i}.f - 2*stats.linReg.gapDur.inspPeakTime{i}.delta, '--', 'color', cmap(i - 1, :)); hold on;
%     label(h, strcat('R^2 = ', num2str(stats.linReg.gapDur.inspPeakTime{i}.r ^2)), 'location', 'west', 'fontsize', 13, 'fontweight', 'bold')
end
%
figure(4)
set(gca, 'tickdir', 'out', 'fontsize', 50)
legend('a-h', 'a-h linear fit', 'a-j', 'a-j linear fit', 'a-i', 'a-i linear fit', 'a-b', 'a-b linear fit')

figure(5)
set(gca, 'tickdir', 'out', 'fontsize', 50)
legend('a-h', 'a-h linear fit', 'a-j', 'a-j linear fit', 'a-i', 'a-i linear fit', 'a-b', 'a-b linear fit')

figure(6)
set(gca, 'tickdir', 'out', 'fontsize', 50)
legend('a-h', 'a-h linear fit', 'a-j', 'a-j linear fit', 'a-i', 'a-i linear fit', 'a-b', 'a-b linear fit')


%% plot the individual trials as different colors given transition type
cmap = colormap('lines');
alpha = 0.5; % to set color of zero line

for i = 2 : length(trans)
    if isempty(breathArray{i}.sylOffset) | length(breathArray{i}.sylOffset(1, :)) == 1
        continue
    end
    
    figure(2)
    x = ([1 : length(breathArray{i}.sylOffset(:, 1))] ./ (fs / 1000)) - (preWindow * 32 / 1000);
    plot(x, breathArray{i}.sylOffset', 'color', cmap(i - 1, :)); hold on

    figure(3)
    h = plot(x, mean(breathArray{i}.sylOffset, 2)); hold on
    X = breathArray{i}.sylOffset;
    ste = std(X') ./ sqrt(length(X(1, :)));
    shadedErrorBar(x', mean(X, 2), ste, 'lineProps', '-', 'lineprops', {'color', cmap(i - 1, :)}); hold on; % standard error bars -- 
%     DO 95% CONFIDENCE INTERVALS INSTEAD -- [sd, confInt] = bootstrap(X)
    label(h, trans(i), 'location', 'east', 'fontsize', 20, 'fontweight', 'bold')
    hold on
    
    if isempty(breathArray{i}.gapOffset)
        continue
    end

    % wrt gap offsets
    figure(4)
    x = ([1 : length(breathArray{i}.gapOffset(:, 1))] ./ (fs / 1000)) - (preWindow * 32 / 1000);
    plot(x, breathArray{i}.gapOffset', 'color', cmap(i - 1, :)); hold on


    figure(5)
    h = plot(x, mean(breathArray{i}.gapOffset, 2)); hold on
    X = breathArray{i}.gapOffset;
    ste = std(X') ./ sqrt(length(X(1, :)));
    shadedErrorBar(x', mean(X, 2), ste, 'lineProps', '-', 'lineprops', {'color', cmap(i - 1, :)}); hold on; % standard error bars -- 
    % DO 95% CONFIDENCE INTERVALS INSTEAD -- [sd, confInt] = bootstrap(X)
    label(h, trans(i), 'location', 'east', 'fontsize', 20, 'fontweight', 'bold')
    hold on
    
    
end

% === wrt gap onset/syllable offset ====

% == song stops! ==
% figure(2)
% x = ([1 : length(breathArray{1}.songStop.air(:, 1))] ./ (fs / 1000)) - (preWindow * 32 / 1000);
% plot(x, breathArray{1}.songStop.air, 'color', cmap(i, :)); hold on

figure(2)
vline(0,'k', 'syllable offset'); hold on;
xlabel('Time (ms)')
ylabel('mV')
% title(strcat('air sac pressure-',trans(1), '_', trans(2), '_', trans(3), '_', trans(4)))
z1 = min(x) : 1 : max(x);
z2 = zeros(length(z1), 1);
plot(z1, z2, 'color',[0,0,0]+alpha)
set(gca, 'tickdir', 'out', 'fontsize', 30)

% figure(3) % plot song stoppings
% h = plot(x, mean(breathArray{1}.songStop.air, 2)); hold on
% X = breathArray{1}.songStop.air;
% ste = std(X') ./ sqrt(length(X(1, :)));
% shadedErrorBar(x', mean(X, 2), ste, 'lineProps', '-', 'lineprops', {'color', cmap(i, :)}); hold on; % standard error bars -- 
% % DO 95% CONFIDENCE INTERVALS INSTEAD -- [sd, confInt] = bootstrap(X)
% label(h, 'song stop', 'location', 'east', 'fontsize', 20, 'fontweight', 'bold')
% hold on

figure(3)
vline(0,'k', 'syllable offset'); hold on;
xlabel('Time (ms)')
ylabel('mV')
% title(strcat('average air sac pressure-', trans(1), '_', trans(2), '_', trans(3), '_', trans(4)))
z1 = min(x) : 1 : max(x);
z2 = zeros(length(z1), 1);
plot(z1, z2, 'color',[0,0,0]+alpha)
set(gca, 'tickdir', 'out', 'fontsize', 30)

% wrt gap offset/onset of subsequent syllable
figure(4)

vline(0,'k', 'gap offset'); hold on;
xlabel('Time (ms)')
ylabel('mV')
% title(strcat('air sac pressure-',trans(1), '_', trans(2), '_', trans(3), '_', trans(4)))
z1 = min(x) : 1 : max(x);
z2 = zeros(length(z1), 1);
plot(z1, z2, 'color',[0,0,0]+alpha)
set(gca, 'tickdir', 'out', 'fontsize', 30)

figure(5)

vline(0,'k', 'gap offset'); hold on;
xlabel('Time (ms)')
ylabel('mV')
% title(strcat('average air sac pressure-', trans(1), '_', trans(2), '_', trans(3), '_', trans(4)))
z1 = min(x) : 1 : max(x);
z2 = zeros(length(z1), 1);
plot(z1, z2, 'color',[0,0,0]+alpha)
set(gca, 'tickdir', 'out', 'fontsize', 30)

%% ==== statistical tests ====
% % do a running statistical test of air sac pressure for the two transitions
% % what's the earliest time that they are statistically different?
% X = breathArray{3}; --- I don't think this is right.... too many points
% seem statistically significant....
% Y = breathArray{4};
% binSize = 1 * (fs / 1000);
% alpha = 0.001;
% [isSig, pVals] = ek_sigHistTime(X, Y, binSize, alpha);
% figure(3)
% scatter((isSig ./ (fs / 1000)) - (window * 32 / 1000), zeros(length(isSig), 1), 'k')
% 
% figure(2)
% scatter((isSig ./ (fs / 1000)) - (window * 32 / 1000), zeros(length(isSig), 1), 'k')

% %% use d' method from Wohlgemuth....
% % d' = (aj - ai) / sqrt((std_aj^2 + std_ai^2) / 2)
% dPrime = (mean(breathArray{3}, 2) - mean(breathArray{4}, 2)) ./ sqrt((std(breathArray{3}').^2 + std(breathArray{4}').^2) ./ 2);
% figure
% plot(breathArray{3}); hold on;
% plot(breathArray{4}); hold on;
% plot(dPrime)

%% use ROC AUC analysis from Lena


%% ==== tidal volumes ====
% inspiratory vs expiratory TVs before and after inspiration

% segement breaths into inspirations and expirations
% below vs above 2.3 V? - test on baseline breathing
% 
% for i = 1 : length(songs)
%     
%     

% can be found in ek_plotChans.m -- calls ek_segmentBreaths.m




%     
%     
%     
% end

%% ======= Convergent syllables! ======= not working properly yet
% Load in data
bird = 'gr98bu75';
date = '081419';
load(strcat(bird, '_airsac_', date, '.mat'))
trans = 'jabh';

% iterate through all songs and find these transitions
% at these transitions put 200 ms after syllable onset into array
for i = 1 : length(trans) % initialize arrays
%     eval(strcat('breathArray.', trans(i), ' = []'));
    breathArray{i} = [];
    breathParams{i}.inspPeakAmp = [];
    breathParams{i}.inspPeakTime = [];
    breathParams{i}.inspVolume = [];
    breathParams{i}.inspDur = [];
    breathParams{i}.gapDur = [];
end

fs = data{1}.syls.Fs;
% figure
for i = 1 : length(data)
    if isempty(data{i})
        continue
    end
    
%     figure
%     subplot(2, 1, 1)
%     plot(data{i}.song); hold on;
%     subplot(2, 1, 2)
%     plot(data{i}.air);
    for j = 2 : length(trans)
    % find transitions
        transInds = strfind(data{i}.syls.labels, strcat(trans(j), trans(1))); % searches for different ordering than above divergent syls
        transSylOnset = data{i}.syls.onsets(transInds + 1);
 
        
    % multiply time point of syllable onset/offset by fs/1000
        transSylOnset = transSylOnset .* (fs / 1000);
    
    % grab 100 ms after syllable offset
        window = 200 * fs/1000;
%         postWindow = 200 * fs/1000;
        
        for k = 1 : length(transInds) - 1
            breathArray{j} = [breathArray{j} data{i}.air(transSylOnset(k) - window : transSylOnset(k) + window)]; % also add in song stopping...
            breathParams{j}.gapDur = [breathParams{j}.gapDur  data{i}.syls.offsets(transInds(k)) - data{i}.syls.onsets(transInds(k) + 1)];
            [Y, I] = min(data{i}.air(transSylOnset(k): transSylOnset(k) + 100));
            breathParams{j}.inspPeakAmp = [breathParams{j}.inspPeakAmp Y];
            breathParams{j}.inspPeakTime = [breathParams{j}.inspPeakTime I];
            breathParams{j}.inspVolume = [breathParams{j}.inspVolume sum(data{i}.air(data{i}.syls.offsets(transInds(k)) : data{i}.syls.onsets(transInds(k) + 1)))];
        end           
    end
end
%% Make figures!
cmap = colormap('lines');
s = 50; % size of scatter points
for i = 2 : length(trans)
    figure(4) % look at relationship between gap duration and peak inspiratory amplitude
    scatter(breathParams{i}.gapDur, breathParams{i}.inspPeakAmp, s, cmap(i - 1, :)); hold on;
    xlabel('Gap Duration (ms)')
    ylabel('Peak Amplitude of Inspiration')
    
    % remove outliers
    [gapDur_noOutliers, outlierIdx] = rmoutliers(breathParams{i}.gapDur);
    inspPeakAmp_noOutliers = breathParams{i}.inspPeakAmp(find(outlierIdx == 0));
    inspVolume_noOutliers = breathParams{i}.inspVolume(find(outlierIdx == 0));
    inspPeakTime_noOutliers = breathParams{i}.inspPeakTime(find(outlierIdx == 0));
    
    % regression to test for correlation
    stats.linReg.gapDur.inspPeakAmp{i}.r = corr(breathParams{i}.gapDur', breathParams{i}.inspPeakAmp');
    [stats.linReg.gapDur.inspPeakAmp{i}.p, stats.linReg.gapDur.inspPeakAmp{i}.S] = polyfit(gapDur_noOutliers,inspPeakAmp_noOutliers,1);
    [stats.linReg.gapDur.inspPeakAmp{i}.f, stats.linReg.gapDur.inspPeakAmp{i}.delta] = polyval(stats.linReg.gapDur.inspPeakAmp{i}.p, gapDur_noOutliers, stats.linReg.gapDur.inspPeakAmp{i}.S');
    h = plot(gapDur_noOutliers, stats.linReg.gapDur.inspPeakAmp{i}.f, 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspPeakAmp{i}.f + (2*stats.linReg.gapDur.inspPeakAmp{i}.delta), '--', 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspPeakAmp{i}.f - (2*stats.linReg.gapDur.inspPeakAmp{i}.delta), '--', 'color', cmap(i - 1, :)); hold on;
    label(h, strcat('rho = ', num2str(stats.linReg.gapDur.inspPeakAmp{i}.r)), 'location', 'east', 'fontsize', 13, 'fontweight', 'bold')
    
    figure(5) % look at relationship between gap duration and inspiratory volume
    scatter(breathParams{i}.gapDur, breathParams{i}.inspVolume, s, cmap(i - 1, :)); hold on;
    xlabel('Gap Duration (ms)')
    ylabel('Inspiratory Volume') % actually do inspiratory volume instead of determining by syl offset/onset
    
     % regression to test for correlation
    stats.linReg.gapDur.inspVolume{i}.r = corr(breathParams{i}.gapDur', breathParams{i}.inspVolume');
    [stats.linReg.gapDur.inspVolume{i}.p, stats.linReg.gapDur.inspVolume{i}.S] = polyfit(gapDur_noOutliers,inspVolume_noOutliers,1);
    [stats.linReg.gapDur.inspVolume{i}.f, stats.linReg.gapDur.inspVolume{i}.delta]  = polyval(stats.linReg.gapDur.inspVolume{i}.p, gapDur_noOutliers, stats.linReg.gapDur.inspVolume{i}.S');
    h = plot(gapDur_noOutliers, stats.linReg.gapDur.inspVolume{i}.f, 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspVolume{i}.f + 2*stats.linReg.gapDur.inspVolume{i}.delta, '--', 'color', cmap(i - 1, :))
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspVolume{i}.f - 2*stats.linReg.gapDur.inspVolume{i}.delta, '--', 'color', cmap(i - 1, :)); hold on;
    label(h, strcat('rho = ', num2str(stats.linReg.gapDur.inspPeakAmp{i}.r)), 'location', 'east', 'fontsize', 13, 'fontweight', 'bold')
    
    figure(6) % look at relationship between gap duration and peak inspiratory time
    scatter(breathParams{i}.gapDur, breathParams{i}.inspPeakTime, s, cmap(i - 1, :)); hold on;
    xlabel('Gap Duration (ms)')
    ylabel('Inspiratory Peak Time From Syllable Offset (ms)')
    
     % regression to test for correlation
    stats.linReg.gapDur.inspPeakTime{i}.r = corr(breathParams{i}.gapDur', breathParams{i}.inspPeakTime');
    [stats.linReg.gapDur.inspPeakTime{i}.p, stats.linReg.gapDur.inspPeakTime{i}.S] = polyfit(gapDur_noOutliers, inspPeakTime_noOutliers,1);
    [stats.linReg.gapDur.inspPeakTime{i}.f, stats.linReg.gapDur.inspPeakTime{i}.delta] = polyval(stats.linReg.gapDur.inspPeakTime{i}.p, gapDur_noOutliers, stats.linReg.gapDur.inspPeakTime{i}.S');
    h = plot(gapDur_noOutliers, stats.linReg.gapDur.inspPeakTime{i}.f, 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspPeakTime{i}.f + 2*stats.linReg.gapDur.inspPeakTime{i}.delta, '--', 'color', cmap(i - 1, :)); hold on;
%     plot(breathParams{i}.gapDur, stats.linReg.gapDur.inspPeakTime{i}.f - 2*stats.linReg.gapDur.inspPeakTime{i}.delta, '--', 'color', cmap(i - 1, :)); hold on;
    label(h, strcat('rho = ', num2str(stats.linReg.gapDur.inspPeakAmp{i}.r)), 'location', 'east', 'fontsize', 13, 'fontweight', 'bold')
end
figure(4)
set(gca, 'tickdir', 'out', 'fontsize', 12)
legend('a-h', 'a-h linear fit', 'a-j', 'a-j linear fit', 'a-i', 'a-i linear fit', 'a-b', 'a-b linear fit')

figure(5)
set(gca, 'tickdir', 'out', 'fontsize', 12)
legend('a-h', 'a-h linear fit', 'a-j', 'a-j linear fit', 'a-i', 'a-i linear fit', 'a-b', 'a-b linear fit')

figure(6)
set(gca, 'tickdir', 'out', 'fontsize', 12)
legend('a-h', 'a-h linear fit', 'a-j', 'a-j linear fit', 'a-i', 'a-i linear fit', 'a-b', 'a-b linear fit')



%% plot the individual trials as different colors given transition type
figure(2)
cmap = colormap('lines');
c = [cmap(2, :); cmap(1, :); cmap(3, :); cmap(4, :)];
alpha = 0.5; % to set color of zero line

for i = 2 : 2 %length(trans)
%     figure(2)
    x = ([1 : length(breathArray{i}(:, 1))] ./ (fs / 1000)) - (window * 32 / 1000);
%     plot(x, breathArray{i}', 'color', c(i, :))
%     hold on
    
    figure(3)
%     transMean.trans(i) = mean(breathArray.(trans(1)).(trans(i)));
    h = plot(x, mean(breathArray{i}, 2)); hold on
    X = breathArray{i};
    ste = std(X') ./ sqrt(length(X(1, :)));
    shadedErrorBar(x', mean(X, 2), ste, 'lineProps', '-'); hold on; % standard error bars -- 
    % DO 95% CONFIDENCE INTERVALS INSTEAD -- [sd, confInt] = bootstrap(X)
    label(h, trans(i), 'location', 'east', 'fontsize', 13, 'fontweight', 'bold')
    hold on
end

figure(2)
vline(0,'k', 'syllable onset'); hold on;
xlabel('Time (ms)')
ylabel('mV')
title(strcat('air sac pressure-',trans(1), '_', trans(2), '_', trans(3), '_', trans(4)))
z1 = min(x) : 1 : max(x);
z2 = zeros(length(z1), 1);
plot(z1, z2, 'color',[0,0,0]+alpha)
set(gca, 'tickdir', 'out', 'fontsize', 30)


figure(3)
vline(0,'k', 'syllable onset'); hold on;
xlabel('Time (ms)')
ylabel('mV')
title(strcat('average air sac pressure-', trans(1), '_', trans(2), '_', trans(3), '_', trans(4)))
z1 = min(x) : 1 : max(x);
z2 = zeros(length(z1), 1);
plot(z1, z2, 'color',[0,0,0]+alpha)
set(gca, 'tickdir', 'out', 'fontsize', 30)









