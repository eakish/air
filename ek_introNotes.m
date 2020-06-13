% find instantaneous respiratory frequency with and without calls
% expiratory and pre/post inspiratory duration
% Load in data

bird = 'or84or24';
date = '120619';
load('or84or24_lRA50uA_airsac_introNotes120619.mat')
fs = data{1, 1}.syls.Fs;
cmap = colormap(lines);

% bird = 'rd8';
% date = '110719';
% load('rd8_airsac_airsac_calls110719.mat')
% fs = data{1, 1}.syls.Fs;
% cmap = colormap(lines);


% % === females ====
% bird = 'or45';
% date = '050819';
% load('or45_airsac_airsac_calls050819.mat')
% fs = data{1, 2}.syls.Fs;
% cmap = colormap(lines);


% === correlation between inspiratory volume/duration/peak and expiratory
% volume - both calls and song? - similar to humans?

%% plot data
for i = 21 : 40 %length(data)
    if isempty(data{1, i})
        continue
    end
    figure
    ax(1) = subplot(2, 1, 1);
%     plot(data{i}.song)
    N = 1024; OVERLAP = 1020; sigma = 3; F_low = 500; F_high = 10000; filter_type = 'hanningffir';
    figs = 1;
    x = 4;
    y=6;
    w=hamming(N);
    SPTH = 300;
    filtsong=pj_bandpass(data{i}.song,fs,500,10000,'butterworth');%'hanningffir'

    [S,F,T] = spectrogram(filtsong,w,OVERLAP,N,fs);
    T = [1 : length(data{i}.air)] ./ (fs / 1000);
    pp = find(abs(S)<=SPTH); %Find entries with very low power + scale up. Makes display nicer.
    S(pp) = SPTH; S = log(abs(S));
    hold on
    imagesc(T,F,S)
    set(gca, 'ylim', [500 10000], 'xlim', [0, max(T)])
    colormap jet

    ax(2) = subplot(2, 1, 2);
    plot(T, data{i}.air, 'k')
    linkaxes(ax, 'x')
end


%% cut out calls w/ window around onset to see single examples
% grab window around call onset
preWin = 1000;
callWin_pre = preWin * fs / 1000;
postWin = 1000;
callWin_post = postWin * fs / 1000;
callThresh = 5000;

inspAmpThresh = -700; expAmpThresh = 500; durThresh = 2000; % for breath segmentation

breathArray_calls = [];
callInd = 0;
for i = 1 : length(data)
    if isempty(data{i})
        continue
    end
    
    % find transitions
    if strcmp(bird, 'rd8')
        callInds = [strfind(data{i}.syls.labels, 'c') strfind(data{i}.syls.labels, 'k') strfind(data{i}.syls.labels, 'i')];
    else
        callInds = strfind(data{i}.syls.labels, 'i');
    end
    callInds = callInds(1); % look at first intro note
    
    callTimes = data{i}.syls.onsets(callInds) .* (fs / 1000); % turn this into call exp onset
    [~, ~, ~,  ~, ~, inspEnd] = ek_segmentBreaths(data{i}.air, bird, durThresh, inspAmpThresh, expAmpThresh);
    c = [];
    for j = 1 : length(callTimes)
        expCall = inspEnd(find(inspEnd < callTimes(j) + 1000));
        c = [c expCall(end)];
    end
    callTimes = c;
%     figure
%     plot(data{i}.air); hold on
%     scatter(c, zeros(length(c), 1), 'r', 'filled')

    % trans times are in samples
    % stim times are in ms
%     if callInd > 50 
%         break
%     end

    for k = 1 : length(callInds)

        if callTimes(k) <= callWin_pre || callTimes(k) + callWin_post > length(data{i}.air)
            continue
        end
        
        
%         if callInd >= 50 
%             break
%         end
        
%         if mod(callInd, 10) == 0
%             figure
%         end

        breaths = data{i}.air(callTimes(k) - callWin_pre : callTimes(k) + callWin_post);
        
        if max(breaths(callWin_pre : callWin_pre + 100 * fs / 1000)) < callThresh % check to make sure there's a call
            continue
        end
        callInd = callInd + 1;
        
        breaths_norm = breaths ./ max(breaths);
        breathArray_calls = [breathArray_calls breaths]; % normalized breaths
        
%         subplot(5, 2, mod(callInd - 1, 10) + 1); hold on
%         ylim([-5000 12000])
%         x = linspace(-preWin, postWin, length(breaths));
%         y = breaths;
%         plot(x, y, 'k', 'linewidth', 1.5)
    end           
end

%% look at trains of calls in relation to breathing

%% compare duration of preceding inspiration with time until next inspiration after call onset - ek 1.9.20
% iterate through call matrix and find breath times
inspAmpThresh = -700; expAmpThresh = 500; durThresh = 2000;
ibi.call= [];
ibi.preCall1 = [];
ibi.preCall2 = [];
ibi.preCall3 = [];
inspDur_preCall1 = [];
callDur = [];
respPeriod_postCall = [];
ibi3.call = [];
ibi3.preCall = [];
callArray = [];
callInd = 0;

% for breath onsets
inspOn.preCall1 = [];
inspOn.postCall1 = [];
expOn.preCall1 = [];
expOn.postCall1 = [];
expOn.call = [];

% for tidal volumes
inspVt.preCall1 = [];
inspVt.postCall1 = [];
expVt.call = [];
expVt.preCall1 = [];
expVt.postCall1 = [];

callBuffer = 1000 * fs / 1000; % filter out calls that have calls occuring in this window around call
callThresh = 5000;

% figure
for i = 1 : length(breathArray_calls(1, :))
    
    % if max amplitude in file is less than expiratory threshold - skip - corrupted data file
    if max(breathArray_calls(i, :)) < expAmpThresh
        continue
    end
    
%     % === if there's a call 500 ms before or after this call, skip. ===
%     % get this from air sac pressure in case I've missed labeling a call...
%     if max(breathArray_calls(callWin_pre - callBuffer : callWin_pre, i)) > callThresh %|| max(breathArray_calls(callWin_pre + 500 : callWin_pre + callBuffer, i)) > callThresh
%         continue % remove breaths from array...
%     end
    
    [breathStartInd, breathEndInd, breathRate,  inspVol, expVol, inspEnd, inspPeak, expPeak, inspDur, expDur] = ek_segmentBreaths(breathArray_calls(:, i), bird, durThresh, inspAmpThresh, expAmpThresh);
    
    if isempty(find(breathStartInd > callWin_pre)) || isempty(find(breathStartInd < callWin_pre))
        continue
    end
    
    % find inspiration preceding call onset
    inspOn_preCall = find(breathStartInd < callWin_pre + 1000);
    inspOn_preCall1 = breathStartInd(inspOn_preCall(end));
%     inspOn_preCall2 = breathStartInd(inspOn_preCall(end - 1));
%     inspOn_preCall3 = breathStartInd(inspOn_preCall(end - 2));
%     inspOn_preCall4 = breathStartInd(inspOn_preCall(end - 3));
%     
%     expOn_preCall = find(inspEnd < callWin_pre + 1000);
%     expOn_preCall1 = inspEnd(expOn_preCall(end));
%     expOn_preCall2 = inspEnd(expOn_preCall(end - 1));
%     expOn_preCall3 = inspEnd(expOn_preCall(end - 2));
%     expOn_preCall4 = inspEnd(expOn_preCall(end - 3));
%     
%     expOn_postCall = find(inspEnd > callWin_pre + 1000);
%     expOn_postCall1 = inspEnd(expOn_postCall(1));
    
%     % find inspiration following call onset
%     inspOn_postCall = find(breathStartInd > callWin_pre + 1000);
%     if length(inspOn_postCall) == 1
%         continue
%     end
%     
%     callArray = [callArray breathArray_calls(:, i)];
%     callInd = callInd + 1;
%     
%     inspOn_postCall1 = breathStartInd(inspOn_postCall(1));
%     inspOn_postCall2 = breathStartInd(inspOn_postCall(2));
%     
%     % find insp & exp vs before call onset... find insp and exp vt (call vt) right after call onset
%     inspVt.preCall1 = [inspVt.preCall1 inspVol(inspOn_preCall(end))];
%     inspVt.postCall1 = [inspVt.postCall1 inspVol(inspOn_postCall(1))];
%     expVt.call = [expVt.call expVol(expOn_preCall(end))];
%     expVt.preCall1 = [expVt.preCall1 expVol(expOn_preCall(end - 1))];
%     expVt.postCall1 = [expVt.postCall1 expVol(expOn_postCall(1))];
%     
%     
%     % insp/exp on times before call
%     inspOn.preCall1 = [inspOn.preCall1 inspOn_preCall1];
%     inspOn.postCall1 = [inspOn.postCall1 inspOn_postCall1];
%     expOn.call = [expOn.call expOn_preCall1];
%     expOn.preCall1 = [expOn.preCall1 expOn_preCall2];
%     expOn.postCall1 = [expOn.postCall1 expOn_postCall1];
%     
%     
%     % find duration of this preceding inspiration
%     inspDur_preCall1 = [inspDur_preCall1 inspDur(inspOn_preCall(end)) * 1000 / fs]; % i don't trust this
%     
%     callDur = [callDur (inspOn_postCall1 - expOn_preCall1) * 1000 / fs];
%     respPeriod_postCall = [respPeriod_postCall (inspOn_postCall2 - inspOn_postCall1) * 1000 / fs];
    
%     % use this to calculate inter-breath interval during call
%     ibi.call = [ibi.call (inspOn_postCall1 - inspOn_preCall1) * 1000 / fs]; % calculate from call exp onset and post call insp onset
%     
%     %  ^ bracketed by breaths
%     ibi3.call = [ibi.call (expOn_postCall1 - expOn_preCall2) * 1000 / fs];
%     ibi3.preCall = [ibi3.preCall (expOn_preCall2 - expOn_preCall4) * 1000 / fs];
%     
% %     % plot these breath trains over each other
% %     if mod(callInd, 10) == 0
% %         figure
% %     end
% %     
% %     subplot(5, 2, mod(callInd - 1, 10) + 1); hold on
% %     T = linspace(0, ibi3.call(end), length(breathArray_calls(expOn_preCall2 : expOn_postCall1, i)));
% %     plot(breathArray_calls(expOn_preCall2 : expOn_postCall1, i), 'r')
% %     
% %     T = linspace(0, ibi3.preCall(end),length(breathArray_calls(expOn_preCall4 : expOn_preCall2, i)));
% %     plot(breathArray_calls(expOn_preCall4 : expOn_preCall2, i), 'k')
% %     scatter(ibi3.preCall(end) * fs / 1000, 0, 'k', 'filled')
% %     scatter(ibi3.call(end) * fs / 1000, 0, 'r', 'filled')
% %     legend('3 expirations with call', '3 expirations before call')
%     
%     % find IBIs of up to 3 inspirations preceding call
%     ibi.preCall1 = [ibi.preCall1 (inspOn_preCall1 - inspOn_preCall2) * 1000 / fs];
%     ibi.preCall2 = [ibi.preCall2 (inspOn_preCall1 - inspOn_preCall3) * 1000 / fs];
%     ibi.preCall3 = [ibi.preCall3 (inspOn_preCall1 - inspOn_preCall4) * 1000 / fs];
end

%% plot stuff
% inspVt.preCall1 = [inspVt.preCall1 inspVol(inspOn_preCall(end))];
%     inspVt.postCall1 = [inspVt.postCall1 inspVol(inspOn_postCall(1))];
%     expVt.call = [expVt.call expVol(expOn_preCall(end))];
%     expVt.preCall1 = [expVt.preCall1 expVol(expOn_preCall(end - 1))];
%     expVt.postCall1 =


figure
scatter(abs(inspVt.preCall1), expVt.call, 'k') % divide up by time of call...
xlabel('Insp_V_t before call')
ylabel('Exp_V_t during call')

figure
scatter(abs(inspVt.preCall1), expVt.preCall1, 'k')
xlabel('Insp_V_t before call')
ylabel('Exp_V_t before call')

figure
scatter(abs(inspVt.postCall1), expVt.call, 'k')
xlabel('Insp_V_t after call')
ylabel('Exp_V_t during call')

%%

figure
scatter(abs(inspVt.preCall1), callDur, 'k')
xlabel('Insp_V_t before call')
ylabel('Call Duration')
set(gca, 'tickdir', 'out', 'fontsize', 15)

%% click on breaths
ek_mouseOverBreaths(abs(inspVt.preCall1), expVt.preCall1, callArray)

%% === breath onsets before and after call ====
[r, pval] = corr(inspDur_preCall1', callDur') % weak positive correlation, not statistically significant -- more data needed
[p, S] = polyfit(inspDur_preCall1,callDur,1);
[f, delta] = polyval(p, inspDur_preCall1, S');
% compare pre-call insp duration to call IBI
figure
scatter(inspDur_preCall1, callDur, 'k')
hold on
h = plot(inspDur_preCall1, f, 'k');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'north', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Insp duration before call (ms)') % doesn't look right.... it's negative and too long
ylabel('Call expiration duration (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

%% compare pre-call exp time to time of following first insp
x = (expOn.preCall1 - callWin_pre) * 1000 / fs;
y = (inspOn.postCall1 - callWin_pre) * 1000 / fs;
[r, pval] = corr(x', y') % weak negative correlation
[p, S] = polyfit(x,y,1);
[f, delta] = polyval(p, x, S');

figure
scatter(x, y, 'k')

hold on
h = plot(x, f, 'k');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'south', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Pre-call exp onset (ms)') % doesn't look right.... it's negative and too long
ylabel('Post-call insp onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

% look only at calls that occur during expirations
expInd = find(x < -200);
x = x(expInd);
y = y(expInd);
[r, pval] = corr(x', y') % weak negative correlation, statistically significant
[p, S] = polyfit(x,y,1);
[f, delta] = polyval(p, x, S');

scatter(x, y, 'r')
h = plot(x, f, 'r');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'west', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Pre-call exp onset (ms)') % doesn't look right.... it's negative and too long
ylabel('Post-call insp onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

legend('call onset during insp', 'call onset during exp?')

%% compare pre-call exp time to time of following first exp
x = (expOn.preCall1 - callWin_pre) * 1000 / fs;
y = (expOn.postCall1 - callWin_pre) * 1000 / fs;

[r, pval] = corr(x', y') % weak negative correlation
[p, S] = polyfit(x,y,1);
[f, delta] = polyval(p, x, S');

figure
scatter(x, y, 'k')

hold on
% h = plot(x, f, 'k');
% label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'south', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Pre-call exp onset (ms)') % doesn't look right.... it's negative and too long
ylabel('Post-call exp onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

% look only at calls that occur during expirations

%
expInd = find(x < -250);
x = x(expInd);
y = y(expInd);

% remove outliers
[x_noOutliers, outlierIdx] = rmoutliers(x);
y_noOutliers = y(find(outlierIdx == 0));

[r, pval] = corr(x_noOutliers', y_noOutliers') % weak negative correlation, statistically significant
[p, S] = polyfit(x_noOutliers,y_noOutliers,1);
[f, delta] = polyval(p, x_noOutliers, S');

scatter(x_noOutliers, y_noOutliers, 'r')
h = plot(x_noOutliers, f, 'r');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'west', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Pre-call exp onset (ms)') % doesn't look right.... it's negative and too long
ylabel('Post-call exp onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

% look only at calls that occur during inspirations
x = (expOn.preCall1 - callWin_pre) * 1000 / fs;
y = (expOn.postCall1 - callWin_pre) * 1000 / fs;

expInd = find(x > -250);
x = x(expInd);
y = y(expInd);

% remove outliers
[x_noOutliers, outlierIdx] = rmoutliers(x);
y_noOutliers = y(find(outlierIdx == 0));

[r, pval] = corr(x_noOutliers', y_noOutliers') % weak negative correlation, statistically significant
[p, S] = polyfit(x_noOutliers,y_noOutliers,1);
[f, delta] = polyval(p, x_noOutliers, S');

scatter(x_noOutliers, y_noOutliers, 'b')
h = plot(x_noOutliers, f, 'b');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'east', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Pre-call exp onset (ms)') % doesn't look right.... it's negative and too long
ylabel('Post-call exp onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)


legend('all calls', 'late call onsets', 'linear fit to early calls no outliers', 'early call onsets', 'linear fit to early calls no outliers')


%% is there a relationship between respiratory volume due to call and subsequent breaths? how long does it take for respiratory system to recover from being perturbed by a call?
figure
scatter(callDur, respPeriod_postCall, 'k')

%% compare inter-breath interval distributions to call IBI - double check that this is right
% figure
% nhist(ibi) % is this real? how does breath segmentation actually look?
% xlabel('Time (ms)')

figure
nhist(ibi3)
title('3 expiration duration')
xlabel('Time (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

%% === call-triggered average ===
x = linspace(-preWin, postWin, length(breaths));

figure
subplot(2, 1, 1)
plot(x, breathArray_calls)
ylim([-5000 12000])
xlabel('Time from call onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

subplot(2, 1, 2)
ste = std(breathArray_calls') ./ sqrt(length(breathArray_calls(1, :)));
shadedErrorBar(x', mean(breathArray_calls, 2), ste, 'lineProps', '-');
ylim([-5000 12000])
xlabel('Time from call onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

%% heatmap of resp aligned to calls sorted by pre-call insp -- normalized to max inspiration to make it easier to look at
[~, inspDur_preCall_sortInd] = sort(inspDur_preCall1);
callArray_sorted = callArray(:, inspDur_preCall_sortInd);
callArray_sorted_norm = (callArray_sorted ./ -min(callArray_sorted));
y = (1 : 1 : length(inspOn.preCall1))';


%     inspOn.preCall1 = [inspOn.preCall1 inspOn_preCall1 - expOn_preCall1];
%     inspOn.postCall1 = [inspOn.postCall1 inspOn_postCall1 - expOn_preCall1];
%     expOn.preCall1 = [expOn.preCall1 expOn_preCall2 - expOn_preCall1];
%     expOn.postCall1 = [expOn.postCall1 expOn_postCall1 - expOn_preCall1];


figure
imagesc(callArray_sorted_norm')
hold on
scatter(inspOn.preCall1(inspDur_preCall_sortInd), y, 'w', 'filled')
scatter(inspOn.postCall1(inspDur_preCall_sortInd), y, 'w', 'filled')
scatter(expOn.preCall1(inspDur_preCall_sortInd), y, 'r', 'filled')
scatter(expOn.postCall1(inspDur_preCall_sortInd), y, 'r', 'filled')

% T = linspace(-preWin,postWin,length(callArray_sorted(:, 1)));
% xticks = 1 : (length(T) - 1) / 6 : length(T);
% xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing normalized to max inspiration')
xlabel('Time relative to call onset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% sort wrt pre-call expiration onset
[~, expOn_preCall_sortInd] = sort(expOn.preCall1);
callArray_sorted = callArray(:, expOn_preCall_sortInd);
callArray_sorted_norm = (callArray_sorted ./ -min(callArray_sorted));
y = (1 : 1 : length(inspOn.preCall1))';


%     inspOn.preCall1 = [inspOn.preCall1 inspOn_preCall1 - expOn_preCall1];
%     inspOn.postCall1 = [inspOn.postCall1 inspOn_postCall1 - expOn_preCall1];
%     expOn.preCall1 = [expOn.preCall1 expOn_preCall2 - expOn_preCall1];
%     expOn.postCall1 = [expOn.postCall1 expOn_postCall1 - expOn_preCall1];


figure
imagesc(callArray_sorted_norm')
hold on
scatter(inspOn.preCall1(expOn_preCall_sortInd), y, 'w', 'filled')
scatter(inspOn.postCall1(expOn_preCall_sortInd), y, 'w', 'filled')
scatter(expOn.preCall1(expOn_preCall_sortInd), y, 'r', 'filled')
scatter(expOn.postCall1(expOn_preCall_sortInd), y, 'r', 'filled')

T = linspace(-preWin,postWin,length(callArray_sorted(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing normalized to max inspiration')
xlabel('Time relative to call onset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar


%% not sorted
callArray_norm = (callArray ./ -min(callArray));
figure
imagesc(callArray_norm')
hold on
y = (1 : 1 : length(inspOn.preCall1))';
scatter(inspOn.preCall1, y, 'w', 'filled')
scatter(inspOn.postCall1, y, 'w', 'filled')
scatter(expOn.preCall1, y, 'r', 'filled')
scatter(expOn.postCall1, y, 'r', 'filled')

T = linspace(-preWin,postWin,length(callArray(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing normalized to max inspiration')
xlabel('Time relative to call onset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar


%% post calls = is post call expiration duration affected by resp phase of call?

% === call triggered average ===
% don't use calls with calls right before



%% breath parameters during calls - look at code below


%% segment breaths, break up into: call breath (pre insp params, exp params, and post insp params)
inspAmpThresh = -500; expAmpThresh = 500; durThresh = 1000;

breathParams.inspPeakAmp = {};
breathParams.expPeakAmp = {};

breathParams.inspVolume = {};
breathParams.expVolume = {};

breathParams.inspDur = {};
breathParams.expDur = {};

breathParams.inspStart = {};
breathParams.inspEnd = {};
breathParams.expEnd = {};
breathParams.isCall = {};

breathInd = 0;
for i = 1 : length(data)
    if isempty(data{i})
        continue
    end 
   
    breathInd = breathInd + 1;
    
    [breathStartInd, breathEndInd, breathRate,  inspVt, expVt, inspEnd, inspPeak, expPeak, inspDur, expDur] = ek_segmentBreaths(data{i}.air, bird, durThresh, inspAmpThresh, expAmpThresh);

    breathParams.inspPeakAmp{breathInd} = inspPeak;
    breathParams.expPeakAmp{breathInd} = expPeak;

    breathParams.inspVolume{breathInd} = inspVt;
    breathParams.expVolume{breathInd} = expVt;
    
    breathParams.inspDur{breathInd} = inspDur ./ (fs * 1000);
    breathParams.expDur{breathInd} = expDur ./ (fs * 1000);
    breathParams.inspStart{breathInd} = breathStartInd;
    breathParams.inspEnd{breathInd} = inspEnd; % find expirations that have calls
    breathParams.expEnd{breathInd} = breathEndInd;
    
    % == see if expiration contains call ==
    % if yes, plot in different color
    callColor = zeros(length(breathStartInd), 3);
    callArray = zeros(length(breathStartInd), 1);
    for j = 1 : length(callArray)
        % iterate through all breaths and see if it contains a call. make array denoting call and no call breaths.
        isCall = intersect(find(data{i}.syls.onsets > breathStartInd), find(data{i}.syls.offsets < breathEndInd));
        if isCall
            callArray(j) = 1;
            callColor(j, :) = cmap(2, :);
        else
            callColor(j, :) = cmap(1, :);
        end
    end
    breathParams.isCall{breathInd} = call;
    
    %  data{i}.syls.onsets data{i}.syls.offsets .labels
    
    % == also look at following inspiration ==
    figure(1)
    scatter(inspPeak, expPeak, callColor)
    hold on

    figure(2)
    scatter(expDur, expPeak, callColor)
    hold on
    
    figure(3)
    scatter(inspDur, expPeak, callColor)
    hold on
    
    figure(4)
    scatter(inspVt, expVt, callColor)
    hold on
    
    figure(5)
    scatter(expPeak, expVt, callColor) % check to make sure that tidal volume is calculated correctly
    hold on
end

figure(1) % incorprate plot breaths
xlabel('inspiratory peak')
ylabel('expiratory peak')

figure(2) % incorprate plot breaths
xlabel('expiratory duration')
ylabel('expiratory peak')

figure(3) % incorprate plot breaths
xlabel('inspiratory duration')
ylabel('expiratory peak')

figure(4) % incorprate plot breaths
xlabel('inspiratory Vt')
ylabel('expiratory Vt')

figure(5) % incorprate plot breaths
xlabel('expiratory peak')
ylabel('expiratory Vt')



%% call rate vs respiratory rate


% === look at how resp phase compares to song onset and how it's effected

for i = 1 : 5
    figure; hold on
    a(1) = subplot(2, 1, 1);
    plot(data{i}.song, 'k')
    a(2) = subplot(2, 1, 2);
    plot(data{i}.air, 'k')
    linkaxes(a, 'x')
end


%% relation of calls to respiratory phase -- how do calculate?
% low-pass filter -- lowpass matlab not working (why?) -- needed to use butterworth filter instead
fc = 25; % cutoff frequency
Wn=2*fc/fs; % the normalized cutoff frequency for butter function
nn=6; % the order
[b,a] = butter(nn,Wn,'low');
breathLow=filter(b,a,data{1}.air);

figure; hold on
% plot(data{1}.air)
plot(breathLow ./ 20, 'r')

%
% four = fft(breathLow);
% n = length(breathLow);
% f = (0:n-1)*(fs/n);
% power = abs(four).^2/n; 
% figure
% plot(f, power)
% xlabel('Frequency')
% ylabel('Power')
%

% hilbert transform to find phase
x = breathLow;
y = hilbert(x);
sigphase = atan2(imag(y),real(y));
% figure
plot(sigphase * 180 / pi, 'k')
% % or
% sigphase = angle(y);
% figure
% plot(sigphase)

% %% == instantaneous frequency === or instfreq(y,fs,'Method','hilbert')
% instfrq = fs/(2*pi)*diff(unwrap(angle(y)));
% figure; hold on
% plot(x ./ 100, 'r')
% plot(instfrq, 'k')


