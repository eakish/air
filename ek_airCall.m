% find instantaneous respiratory frequency with and without calls
% expiratory and pre/post inspiratory duration
% Load in data
% 
bird = 'or84or24';
date = '121819';
% load('or84or24_lRA50uA_airsac_stim_calls120619.mat')
load('or84or24_lRA50uA_airsac_stim_calls_noStim120619.mat')
% load('or84or24_lRA50uA_airsac_calls120619.mat')
% load('or84or24_lRA90uA_airsac_stim_noCall121819.mat')
fs = 32000; %data{1, 1}.syls.Fs;
cmap = colormap(lines);
baselineSub = 0;
callThresh = 5000;
maxAmp = 15000;

% bird = 'rd8';
% date = '110719';
% load('rd8_airsac_airsac_calls110719.mat')
% fs = 32000; %data{1, 1}.syls.Fs;
% cmap = colormap(lines);
% baselineSub = -800;
% callThresh = 7000;
% maxAmp = 15000;

% bird = 'gr98bu75';
% date = '110719';
% cd '/home/eszter/Data/airsac/males/gr98bu75/preProcessed/';
% load('gr98bu75_airsac_airsac_calls081419.mat')
% fs = 32000; %data{1, 1}.syls.Fs;
% cmap = colormap(lines);
% baselineSub = -1000; % how much to subtract from traces to get a baseline of zero --- this is currently just based on looking at traces and is pretty arbitrary. the capacitor filtering is trash.
% callThresh = 4000;

% % === females ====
% bird = 'or45';
% date = '050819';
% load('or45_airsac_airsac_calls050819.mat')
% fs = data{1, 2}.syls.Fs;
% cmap = colormap(lines);


% === correlation between inspiratory volume/duration/peak and expiratory
% volume - both calls and song? - similar to humans?

%% plot data
for i = 1 :20 %length(data)
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
    plot(T, data{i}.air - baselineSub, 'k')
    linkaxes(ax, 'x')
end


%% cut out calls w/ window around onset to see single examples
% grab window around call onset
preWin = 2000;
callWin_pre = preWin * fs / 1000;
postWin = 1000;
callWin_post = postWin * fs / 1000;
% callThresh = 5000;

inspAmpThresh = -1000; expAmpThresh = 500; durThresh = 2000; % for breath segmentation

breathArray_calls = [];
callInd = 0;

% breathParams.noCalls = {};

breathParams.noCalls.vol.insp = [];
breathParams.noCalls.vol.exp = [];
breathParams.noCalls.inspEnd  = [];
breathParams.noCalls.peak.insp = [];
breathParams.noCalls.peak.exp = [];
breathParams.noCalls.dur.insp = [];
breathParams.noCalls.dur.exp = [];
breathParams.noCalls.phase.timeFromPrevExp = [];
breathParams.noCalls.phase.timeFromPrevInsp = [];
breathParams.noCalls.phase.timeBetweenExpAndNextInsp = [];

callParams.interCallinterval.single = [];
callParams.interCallinterval.all = [];
callSeq.autoCorr.air = {};
callSeq.autoCorr.sound = {};
callSeq.crossCorr.soundVSair = {};
callSeq.crossCorr.lags = {};
callSeq.autoCorr.lags = {};


% == for hamish/phase analysis ==
forHamish(1).breathing = [];
forHamish(1).callOnset = [];
hamish = 0;

noCallBreathing = [];
for i = 1 : length(data)
    if isempty(data{i})
        continue
    end
    
%     % === auto-correlation of breathing and sound envelopes to find rhythmic structure ===
%     [air_corr, lags]  = xcorr(data{i}.air);
%     callSeq.autoCorr.air{i} = air_corr;
%     callSeq.autoCorr.lags{i} = lags;
%     
%      % use PJ method to find amplitude envelopes
%     N = 1024; OVERLAP = 1020; sigma = 3; F_low = 500; F_high = 10000; filter_type = 'hanningffir';
%     filtSong = highpass(data{i}.song, F_low, fs);
%     filtSong = lowpass(filtSong, F_high, fs);
%     filtSong = rms(filtSong, 2);
%     win =  2.5 * fs / 1000;
%     envSong = smooth(filtSong, win);
    
% 
%     callSeq.autoCorr.sound{i} = xcorr(envSong);
%     
%     % == cross-correlation of normalized sound vs airsac pressure ==
%     [callSeq.crossCorr.soundVSair{i}, callSeq.crossCorr.lags{i}] = xcorr(envSong ./ max(envSong), data{i}.air ./ max(data{i}.air));
    
    
    % === find calls ===
    air = data{i}.air - baselineSub;
    [breathStartInd, breathEndInd, ~, ~, ~, inspEnd, inspPeak, expPeak, inspDur, expDur] = ek_segmentBreaths(air, bird, durThresh, inspAmpThresh, expAmpThresh);
    c = [];
    for j = 1 : length(inspEnd)
        if (expPeak(j) > callThresh) && (expDur(j) > 100 * fs / 1000)
            callExpDur = length([find(air(breathStartInd : breathEndInd) > 2000)]);
            if callExpDur < 1000 %%&& expPeak(j) > callThresh && callExpDur > 500
                c = [c inspEnd(j)];
            end
        end
    end
    callTimes = c;
    
%     if ~isempty(c)
%             figure; plot(envSong); hold on
%             plot(air)
%             x = 1;
%     end
    
    
    % === get rid of bad quality recording files ===
    if min(air(c - 3000 : c - 2000) > -500)
        continue
    end
    
    if max(air) > maxAmp
        continue
    end
    
    
    % === grab single call file ===
    if length(callTimes) == 1 %&& isempty(forHamish.breathing)
        hamish = hamish + 1;
        forHamish(hamish).callOnset = callTimes;
        forHamish(hamish).breathing = data{i}.air;
    end
    
    for j = 1 : length(callTimes)
        
        if j > 1
            callParams.interCallinterval.single = [callParams.interCallinterval.single callTimes(j) - callTimes(j - 1)];
        end
        
        % === find the difference between this call and all the other ones ===
        for k = 1 : length(callTimes)
            callParams.interCallinterval.all = [callParams.interCallinterval.all callTimes(j) - callTimes(k)];
        end 
    end
    
%     figure
%     plot(data{i}.air); hold on
%     scatter(c, zeros(length(c), 1), 'r', 'filled')

    % trans times are in samples
    % stim times are in ms
%     if callInd > 50 
%         break
%     end

    for k = 1 : length(callTimes)

        if callTimes(k) <= callWin_pre || callTimes(k) + callWin_post > length(air)
            continue
        end
        
        
%         if callInd >= 50 
%             break
%         end
        
%         if mod(callInd, 10) == 0
%             figure
%         end

        breaths = air(callTimes(k) - callWin_pre : callTimes(k) + callWin_post);
        
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
    
    % === find time periods with no large deflections in breaths right before/after ===
    [breathStartInd, breathEndInd, breathRate,  inspVol, expVol, inspEnd, inspPeak, expPeak, inspDur, expDur] = ek_segmentBreaths(air, bird, durThresh, inspAmpThresh, expAmpThresh);
    
    % find breath params in these no call time periods
%     breathParams.noCalls{i}.inspVol = [];
%     breathParams.noCalls{i}.expVol = [];
%     breathParams.noCalls{i}.inspEnd  = [];
%     breathParams.noCalls{i}.inspPeak = [];
%     breathParams.noCalls{i}.expPeak = [];
%     breathParams.noCalls{i}.inspDur = [];
%     breathParams.noCalls{i}.expDur = [];
%     breathParams.noCalls{i}.breath = [];
    for k = 2 : length(breathStartInd) - 2
        if max(air(breathStartInd(k - 1) :  breathEndInd(k + 1))) > callThresh
            continue
        end
        breathParams.noCalls.vol.insp = [breathParams.noCalls.vol.insp inspVol(k)];
        breathParams.noCalls.vol.exp = [breathParams.noCalls.vol.exp expVol(k)];
        breathParams.noCalls.inspEnd  = [breathParams.noCalls.inspEnd inspEnd(k)];
        breathParams.noCalls.peak.insp = [breathParams.noCalls.peak.insp inspPeak(k)];
        breathParams.noCalls.peak.exp = [breathParams.noCalls.peak.exp expPeak(k)];
        breathParams.noCalls.dur.insp = [breathParams.noCalls.dur.insp inspDur(k) * 1000 / fs];
        breathParams.noCalls.dur.exp = [breathParams.noCalls.dur.exp expDur(k) * 1000 / fs];
        breathParams.noCalls.phase.timeFromPrevExp = [breathParams.noCalls.phase.timeFromPrevExp (inspEnd(k) - inspEnd(k - 1)) * 1000 / fs];
        breathParams.noCalls.phase.timeBetweenExpAndNextInsp = [breathParams.noCalls.phase.timeBetweenExpAndNextInsp (breathStartInd(k + 1) - inspEnd(k - 1)) * 1000 / fs];
        breathParams.noCalls.phase.timeFromPrevInsp = [breathParams.noCalls.phase.timeFromPrevInsp (breathStartInd(k) - breathStartInd(k - 1)) * 1000 / fs];
        
%         breathParams{i}.noCalls.inspVol = [breathParams.noCalls{i}.inspVol inspVol(k)];
%         breathParams{i}.noCalls.expVol = [breathParams.noCalls{i}.expVol expVol(k)];
%         breathParams{i}.noCalls.inspEnd  = [breathParams.noCalls{i}.inspEnd inspEnd(k)];
%         breathParams{i}.noCalls.inspPeak = [breathParams.noCalls{i}.inspPeak inspPeak(k)];
%         breathParams{i}.noCalls.expPeak = [breathParams.noCalls{i}.expPeak expPeak(k)];
%         breathParams{i}.noCalls.inspDur = [breathParams.noCalls{i}.inspDur inspDur(k)];
%         breathParams{i}.noCalls.expDur = [breathParams.noCalls{i}.expDur expDur(k)];
%         breathParams.noCalls{i}.breath = [breathParams.noCalls{i}.breath data{i}.air(breathStartInd(k) : breathEndInd(k))'];
    end
    
end

% save stuff here
save(strcat(bird, '_', date, '_airCallsArray.mat'), 'preWin', 'postWin', 'callWin_pre', 'breathArray_calls', 'breathParams', 'fs', '-v7.3'); % insert variables to save

%% compare duration of preceding inspiration with time until next inspiration after call onset - ek 1.9.20
% iterate through call matrix and find breath times
% load('rd8_110719_airCallsArray.mat')
% fn = 'or84or24_121819_airCallsArray.mat';
fn = 'rd8_110719_airCallsArray.mat';

load(fn)
fileSplit = strsplit(fn, '_');
bird = fileSplit{1};
date = fileSplit{2};

inspAmpThresh = -700; expAmpThresh = 500; durThresh = 2000;
ibi.call= [];
ibi.preCall1 = [];
ibi.preCall2 = [];
ibi.preCall3 = [];
inspDura.preCall1 = [];
expDura.preCall1 = [];
callDur = [];
respPeriod_postCall = [];
ibi3.call = [];
ibi3.preCall = [];
callArray = [];
callInd = 0;

% for breath onsets
inspOn.preCall1 = [];
inspDura.preCall2 = [];
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

callParams.peakAmp = [];

callBuffer = 1000 * fs / 1000; % filter out calls that have calls occuring in this window around call
% callThresh = 5000;

% figure
for i = 1 : length(breathArray_calls(1, :))
    
    % if max amplitude in file is less than expiratory threshold - skip - corrupted data file
    if max(breathArray_calls(i, :)) < expAmpThresh % at baseline... w/o calls
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
    inspOn_preCall2 = breathStartInd(inspOn_preCall(end - 1));
    inspOn_preCall3 = breathStartInd(inspOn_preCall(end - 2));
%     inspOn_preCall4 = breathStartInd(inspOn_preCall(end - 3));
    
    expOn_preCall = find(inspEnd < callWin_pre + 1000);
    expOn_preCall1 = inspEnd(expOn_preCall(end));
    expOn_preCall2 = inspEnd(expOn_preCall(end - 1));
    expOn_preCall3 = inspEnd(expOn_preCall(end - 2));
%     expOn_preCall4 = inspEnd(expOn_preCall(end - 3));
    
    expOn_postCall = find(inspEnd > callWin_pre + 1000);
    expOn_postCall1 = inspEnd(expOn_postCall(1));
    
    % find inspiration following call onset
    inspOn_postCall = find(breathStartInd > callWin_pre + 1000);
    if length(inspOn_postCall) == 1
        continue
    end
    
    callArray = [callArray breathArray_calls(:, i)];
    callInd = callInd + 1;
    
    inspOn_postCall1 = breathStartInd(inspOn_postCall(1));
    inspOn_postCall2 = breathStartInd(inspOn_postCall(2));
    
    % find insp & exp vs before call onset... find insp and exp vt (call vt) right after call onset
    inspVt.preCall1 = [inspVt.preCall1 inspVol(inspOn_preCall(end))];
    inspVt.postCall1 = [inspVt.postCall1 inspVol(inspOn_postCall(1))];
    expVt.call = [expVt.call expVol(expOn_preCall(end))];
    expVt.preCall1 = [expVt.preCall1 expVol(expOn_preCall(end - 1))];
    expVt.postCall1 = [expVt.postCall1 expVol(expOn_postCall(1))];
    
    
    % insp/exp on times before call
    inspOn.preCall1 = [inspOn.preCall1 inspOn_preCall1];
    inspOn.postCall1 = [inspOn.postCall1 inspOn_postCall1];
    expOn.call = [expOn.call expOn_preCall1];
    expOn.preCall1 = [expOn.preCall1 expOn_preCall2];
    expOn.postCall1 = [expOn.postCall1 expOn_postCall1];
    
    
    % find duration of this preceding inspiration
    inspDura.preCall1 = [inspDura.preCall1 inspDur(inspOn_preCall(end)) * 1000 / fs];
    inspDura.preCall2 = [inspDura.preCall2 inspDur(inspOn_preCall(end - 1)) * 1000 / fs]; % corresponds to inspiration before the expiration below
    expDura.preCall1 = [expDura.preCall1 expDur(expOn_preCall(end - 1)) * 1000 / fs]; % the last expiration is the call, end - 1 is the expiration preceding the call
    
    callDur = [callDur (inspOn_postCall1 - expOn_preCall1) * 1000 / fs];
    callParams.peakAmp = [callParams.peakAmp max(breathArray_calls(expOn_preCall1 : inspOn_postCall1, i))];
    % === find amplitudes of the two peaks? ===
    
    respPeriod_postCall = [respPeriod_postCall (inspOn_postCall2 - inspOn_postCall1) * 1000 / fs];
    
%     % use this to calculate inter-breath interval during call
%     ibi.call = [ibi.call (inspOn_postCall1 - inspOn_preCall1) * 1000 / fs]; % calculate from call exp onset and post call insp onset
%     
%     %  ^ bracketed by breaths
%     ibi3.call = [ibi.call (expOn_postCall1 - expOn_preCall2) * 1000 / fs];
%     ibi3.preCall = [ibi3.preCall (expOn_preCall2 - expOn_preCall4) * 1000 / fs];
    
%     % plot these breath trains over each other
%     if mod(callInd, 10) == 0
%         figure
%     end
%     
%     subplot(5, 2, mod(callInd - 1, 10) + 1); hold on
%     T = linspace(0, ibi3.call(end), length(breathArray_calls(expOn_preCall2 : expOn_postCall1, i)));
%     plot(breathArray_calls(expOn_preCall2 : expOn_postCall1, i), 'r')
%     
%     T = linspace(0, ibi3.preCall(end),length(breathArray_calls(expOn_preCall4 : expOn_preCall2, i)));
%     plot(breathArray_calls(expOn_preCall4 : expOn_preCall2, i), 'k')
%     scatter(ibi3.preCall(end) * fs / 1000, 0, 'k', 'filled')
%     scatter(ibi3.call(end) * fs / 1000, 0, 'r', 'filled')
%     legend('3 expirations with call', '3 expirations before call')
    
%     % find IBIs of up to 3 inspirations preceding call
%     ibi.preCall1 = [ibi.preCall1 (inspOn_preCall1 - inspOn_preCall2) * 1000 / fs];
%     ibi.preCall2 = [ibi.preCall2 (inspOn_preCall1 - inspOn_preCall3) * 1000 / fs];
%     ibi.preCall3 = [ibi.preCall3 (inspOn_preCall1 - inspOn_preCall4) * 1000 / fs];
end

% save pre-processed data
% cd('preProcessed2')

save(strcat(bird, '_', date, '_airCallsParams'), 'ibi', 'inspDura', 'expDura', 'callDur', 'respPeriod_postCall', 'ibi3', 'callArray', 'inspOn', 'expOn', 'inspVt', 'expVt', 'callParams', '-V7.3');

%% plot stuff
clear all

% bird = 'or84or24';
% date = '121819';

bird = 'rd8';
date = '110719';

load(strcat(bird, '_', date, '_airCallsParams'))
load(strcat(bird, '_', date, '_airCallsArray.mat'))

% === load pre-processed data here ===

% inspVt.preCall1 = [inspVt.preCall1 inspVol(inspOn_preCall(end))];
%     inspVt.postCall1 = [inspVt.postCall1 inspVol(inspOn_postCall(1))];
%     expVt.call = [expVt.call expVol(expOn_preCall(end))];
%     expVt.preCall1 = [expVt.preCall1 expVol(expOn_preCall(end - 1))];
%     expVt.postCall1 =


figure
subplot(3, 1, 1)
scatter(abs(inspVt.preCall1), expVt.call, 'k') % divide up by time of call...
xlabel('Insp_V_t before call')
ylabel('Exp_V_t during call')
set(gca, 'tickdir', 'out', 'fontsize', 15)

subplot(3, 1, 2)
scatter(abs(inspVt.preCall1), expVt.preCall1, 'k')
xlabel('Insp_V_t before call')
ylabel('Exp_V_t before call')
set(gca, 'tickdir', 'out', 'fontsize', 15)

subplot(3, 1, 3)
scatter(abs(inspVt.postCall1), expVt.call, 'k')
xlabel('Insp_V_t after call')
ylabel('Exp_V_t during call')
set(gca, 'tickdir', 'out', 'fontsize', 15)

%% absolute Vt pre call vs duration of call
absVt = inspVt.preCall1 + expVt.preCall1;

[x_noOutliers, outlierIdx] = rmoutliers(absVt);
y_noOutliers = callDur(find(outlierIdx == 0));

[r, pval] = corr(x_noOutliers', y_noOutliers') % weak positive correlation, not statistically significant -- more data needed
[p, S] = polyfit(x_noOutliers, y_noOutliers, 1);
[f, delta] = polyval(p, x_noOutliers, S');

figure
scatter(x_noOutliers, y_noOutliers, 'k'); hold on
h = plot(x_noOutliers, f, 'k');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'north', 'fontsize', 11, 'fontweight', 'bold')

xlabel('Absolute Vt before call')
ylabel('Call Duration (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

ek_mouseOverBreaths(x_noOutliers, y_noOutliers, callArray)

%% absolute Vt pre call vs. call Vt
absVt = inspVt.preCall1 + expVt.preCall1;

[x_noOutliers, outlierIdx] = rmoutliers(absVt);
y_noOutliers = expVt.call(find(outlierIdx == 0));

[r, pval] = corr(x_noOutliers', y_noOutliers') % weak positive correlation, not statistically significant -- more data needed
[p, S] = polyfit(x_noOutliers, y_noOutliers, 1);
[f, delta] = polyval(p, x_noOutliers, S');

figure
scatter(x_noOutliers, y_noOutliers, 'k'); hold on
h = plot(x_noOutliers, f, 'k');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'north', 'fontsize', 11, 'fontweight', 'bold')

xlabel('Absolute Vt before call')
ylabel('Call Vt')
set(gca, 'tickdir', 'out', 'fontsize', 15)

ek_mouseOverBreaths(x_noOutliers, y_noOutliers, callArray)

%% call dur vs. call Vt

[x_noOutliers, outlierIdx] = rmoutliers(callDur);
y_noOutliers = expVt.call(find(outlierIdx == 0));

[r, pval] = corr(x_noOutliers', y_noOutliers') % weak positive correlation, not statistically significant -- more data needed
[p, S] = polyfit(x_noOutliers, y_noOutliers, 1);
[f, delta] = polyval(p, x_noOutliers, S');

figure
scatter(x_noOutliers, y_noOutliers, 'k'); hold on
h = plot(x_noOutliers, f, 'k');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'north', 'fontsize', 11, 'fontweight', 'bold')

xlabel('Call Duration')
ylabel('Call Vt')
set(gca, 'tickdir', 'out', 'fontsize', 15)

ek_mouseOverBreaths(x_noOutliers, y_noOutliers, callArray)

%% call peak amplitude
figure
scatter(expVt.call, callParams.peakAmp, 'k')
xlabel('Exp_V_t Call')
ylabel('Call Peak Amplitude')
set(gca, 'tickdir', 'out', 'fontsize', 15)

figure
scatter(expVt.call, callDur, 'k')
xlabel('Exp_V_t Call')
ylabel('Call Duration')
set(gca, 'tickdir', 'out', 'fontsize', 15)

figure
scatter(callParams.peakAmp, callDur, 'k')
xlabel('Call Peak Amplitude')
ylabel('Call Duration')
set(gca, 'tickdir', 'out', 'fontsize', 15)

figure
scatter3(callDur, expVt.call, callParams.peakAmp, 'k')
xlabel('Call Duration')
ylabel('Exp_V_t Call')
zlabel('Call Peak Amplitude')
set(gca, 'tickdir', 'out', 'fontsize', 15)


%%
figure
scatter3(callDur, abs(inspVt.preCall1), callParams.peakAmp, 'k')
xlabel('Call Duration')
ylabel('Pre-call Insp Vt')
zlabel('Call Peak Amplitude')
set(gca, 'tickdir', 'out', 'fontsize', 15)

%% call dur vs. pre call insp and exp Vt
figure
scatter3(abs(inspVt.preCall1), expVt.preCall1, callDur, 'k')
xlabel('Insp_V_t before call')
ylabel('Exp_V_t before')
zlabel('Call Duration (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

%% click on breaths
ek_mouseOverBreaths(absVt, callDur, callArray)

%% === breath onsets before and after call ====
[r, pval] = corr(inspDura.preCall1', callDur') % weak positive correlation, not statistically significant -- more data needed
[p, S] = polyfit(inspDura.preCall1,callDur,1);
[f, delta] = polyval(p, inspDura.preCall1, S');
% compare pre-call insp duration to call IBI
figure
scatter(inspDura.preCall1, callDur, 'k')
hold on
h = plot(inspDura.preCall1, f, 'k');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'north', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Insp duration before call (ms)') % doesn't look right.... it's negative and too long
ylabel('Call expiration duration (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

%% === exp/insp onset times relative to call ===
breathTime.expPreCall = (expOn.preCall1 - callWin_pre) * 1000 / fs;
breathTime.expPostCall = (expOn.postCall1 - callWin_pre) * 1000 / fs;
breathTime.inspPreCall = (inspOn.preCall1 - callWin_pre) * 1000 / fs;
breathTime.inspPostCall = (inspOn.postCall1 - callWin_pre) * 1000 / fs;
figure; nhist(breathTime)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Time relative to call onset (ms)')



%% === expiratory phase at baseline vs with call ===
expPhase.baseline = breathParams.noCalls.phase.timeFromPrevExp;
expPhase.call = abs((expOn.preCall1 - callWin_pre) * 1000 / fs); % phase of breathing based on exp dur when call occurs

figure; nhist(expPhase)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Time between expiration onsets (ms)')

%% === insp duration vs following expiratory phase === 
clear expPhase
expPhase.baseline = breathParams.noCalls.phase.timeFromPrevExp;
preInspDura.baseline = breathParams.noCalls.dur.insp; % happens right before these expirations --- check if correct ---
expPhase_inspOn.baseline = breathParams.noCalls.phase.timeBetweenExpAndNextInsp;

expPhase.call = (expOn.postCall1 - expOn.preCall1) * 1000 / fs; % phase of breathing wrt to expiratory onset encompassing call
expPhase_inspOn.call = (inspOn.postCall1 - expOn.preCall1) * 1000 / fs;
preInspDura.call = inspDura.preCall2; % pre call 2

figure
subplot(2, 1, 1)
scatter(preInspDura.baseline(1 : end - 1), expPhase.baseline(2 : end), 'k'); hold on
scatter(preInspDura.call, expPhase.call, 'r')
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Duration of Inspiration preceding expiration(ms)')
ylabel('Time between two expiratory onsets (ms)')
legend('Baseline breathing', 'Call separates the expirations')

subplot(2, 1, 2)
nhist(expPhase)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Time between two expiratory onsets (ms)')


%% === with side histograms ===

clear group

preInspDura.baseline = breathParams.noCalls.dur.insp; 
expPhase.baseline =breathParams.noCalls.phase.timeFromPrevExp;
expPhase_inspOn.baseline = breathParams.noCalls.phase.timeBetweenExpAndNextInsp;

preInspDura.call = inspDura.preCall2; % pre call 2
expPhase.call = (expOn.postCall1 - expOn.preCall1) * 1000 / fs; % phase of breathing wrt to expiratory onset encompassing call
expPhase_inspOn.call = (inspOn.postCall1 - expOn.preCall1) * 1000 / fs;

x = [preInspDura.baseline(1 : end - 1) preInspDura.call];
y = [expPhase.baseline(2 : end) expPhase.call];
for i = 1 : length(preInspDura.baseline(1 : end - 1))
    group{i, 1} = 'Baseline';
end
for i = length(preInspDura.baseline(1 : end - 1)) + 1 : length(preInspDura.baseline(1 : end - 1))+ length(preInspDura.call)
    group{i, 1} = 'Call';
end

figure % 'Kernel', 'on', 
scatterhist(x, y, 'Group', group,'Kernel', 'on','Location', 'SouthEast', 'Direction', 'out', 'Color', 'kr', 'LineWidth', [3, 3, 3]);
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Duration of Inspiration preceding expiration(ms)')
ylabel('Time between two expiratory onsets (ms)')
legend('Baseline breathing', 'Separated by call')





%% === insp duration vs following expiratory phase relative to insp onset ===
figure
subplot(2, 1, 1)
scatter(preInspDura.baseline(1 : end - 1), expPhase_inspOn.baseline(2 : end), 'k'); hold on
scatter(preInspDura.call, expPhase_inspOn.call, 'r')
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Duration of Inspiration preceding expiration(ms)')
ylabel('Time between exp onsets and onset of next inspiration(ms)')
legend('Baseline breathing', 'Separated by call')

subplot(2, 1, 2)
nhist(expPhase_inspOn)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Time between two expiratory onsets (ms)')

%% === with side histograms ===
clear group
x = [preInspDura.baseline(1 : end - 1) preInspDura.call];
y = [expPhase_inspOn.baseline(2 : end) expPhase_inspOn.call];
for i = 1 : length(preInspDura.baseline(1 : end - 1))
    group{i, 1} = 'Baseline';
end
for i = length(preInspDura.baseline(1 : end - 1)) + 1 : length(preInspDura.baseline(1 : end - 1))+ length(preInspDura.call)
    group{i, 1} = 'Call';
end

figure % 'Kernel', 'on', 
scatterhist(x, y, 'Group', group, 'Kernel', 'on', 'Location', 'SouthEast', 'Direction', 'out', 'Color', 'kr', 'LineWidth', [3, 3, 3]);
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Duration of Inspiration preceding expiration (ms)')
ylabel('Time between exp onsets and onset of next inspiration (ms)')
legend('Baseline breathing', 'Separated by call')

%% === predictive linear regression for phase shift ===
x_base = preInspDura.baseline(1 : end - 1);
y_base = expPhase_inspOn.baseline(2 : end);

x_call = preInspDura.call;
y_call = expPhase_inspOn.call;

[p, S] = polyfit(x_base, y_base,1);
[f, delta] = polyval(p, x_call, S');

unity = [0 : max(y_call)];
figure; scatter(y_call, f, 'k'); hold on
% plot(unity, unity, 'r')

%% pre call exp onset vs call amplitude
figure
scatter(callParams.peakAmp, breathTime.expPreCall, 'k')
xlabel('Pre-Call Expiratory Onset')
ylabel('Call Peak Amplitude')
set(gca, 'tickdir', 'out', 'fontsize', 15)
hold on
[r, pval] = corr(callParams.peakAmp',  breathTime.expPreCall') % weak negative correlation, statistically significant
[p, S] = polyfit(callParams.peakAmp, breathTime.expPreCall,1);
[f, delta] = polyval(p, callParams.peakAmp, S');
h = plot(callParams.peakAmp, f, 'r');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'west', 'fontsize', 11, 'fontweight', 'bold')

%% === baseline respiratory parameters ===
% breathParams.noCalls.vol.insp = [];
% breathParams.noCalls.vol.exp = [];
% breathParams.noCalls.inspEnd  = [];
% breathParams.noCalls.peak.insp = [];
% breathParams.noCalls.peak.exp = [];
% breathParams.noCalls.dur.insp = [];
% breathParams.noCalls.dur.exp = [];

figure
subplot(3, 1, 1)
nhist(breathParams.noCalls.dur)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Duration (ms)')
title('No call breath durations')

% figure
subplot(3, 1, 2)
nhist(breathParams.noCalls.vol)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Volume (A.U.)')
title('No call breath volumes')

% figure
subplot(3, 1, 3)
nhist(breathParams.noCalls.peak)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Peak amplitude (A.U.)')
title('No call breath peak')

%% === baseline insp/exp duration vs. pre-call insp/exp duration ===
clear expDuration inspDuration

inspDuration.noCall = breathParams.noCalls.dur.insp;
inspDuration.preCall = inspDura.preCall1;
inspDuration.postCall = (expOn.postCall1 - inspOn.postCall1) * 1000 / fs;

figure
% subplot(2, 1, 1)
nhist(inspDuration)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Time (ms)')
title('Inspiratory Durations')
%%
clear expDuration

expDuration.expDur_noCall = breathParams.noCalls.dur.exp;
expDuration.expDur_call = callDur;
expDuration.expDur_preCall = expDura.preCall1;

figure
nhist(expDuration)
set(gca, 'tickdir', 'out', 'fontsize', 15)
xlabel('Time (ms)')
title('Baseline Expiratory Durations vs. Expiratory onset time to call')


%% === baseline insp/exp parameters ===
figure
scatter(breathParams.noCalls.vol.insp, breathParams.noCalls.vol.exp, 'k')
xlabel('Inspiratory volume - leading')
ylabel('Expiratory volume - following')

figure
scatter(breathParams.noCalls.vol.insp(2 : end), breathParams.noCalls.vol.exp(1 : end - 1), 'k')
xlabel('Inspiratory volume - following')
ylabel('Expiratory volume - leading')

figure
scatter(breathParams.noCalls.dur.insp(2 : end), breathParams.noCalls.dur.exp(1 : end - 1), 'k')
xlabel('Inspiratory duration - following') % which is more strongly correlated?
ylabel('Expiratory duration - leading')
hold on
[r, pval] = corr(breathParams.noCalls.dur.insp(2 : end)', breathParams.noCalls.dur.exp(1 : end - 1)') % weak negative correlation, statistically significant
[p, S] = polyfit(breathParams.noCalls.dur.insp(2 : end),breathParams.noCalls.dur.exp(1 : end - 1),1);
[f, delta] = polyval(p, breathParams.noCalls.dur.insp(2 : end), S');
h = plot(breathParams.noCalls.dur.insp(2 : end), f, 'r');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'west', 'fontsize', 11, 'fontweight', 'bold')

figure
scatter(breathParams.noCalls.dur.insp, breathParams.noCalls.dur.exp, 'k')
xlabel('Inspiratory duration - leading')
ylabel('Expiratory duration - following')
hold on
[r, pval] = corr(breathParams.noCalls.dur.insp', breathParams.noCalls.dur.exp') % weak negative correlation, statistically significant
[p, S] = polyfit(breathParams.noCalls.dur.insp,breathParams.noCalls.dur.exp,1);
[f, delta] = polyval(p, breathParams.noCalls.dur.insp, S');
h = plot(breathParams.noCalls.dur.insp, f, 'r');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'west', 'fontsize', 11, 'fontweight', 'bold')

%% === compare pre-call exp time to time of following first insp ===
x = (expOn.preCall1 - callWin_pre) * 1000 / fs;
y = (inspOn.postCall1 - callWin_pre) * 1000 / fs;
[r, pval] = corr(x', y') % weak negative correlation
[p, S] = polyfit(x,y,1);
[f, delta] = polyval(p, x, S');

figure
scatter(abs(x), y, 'k')

hold on
% h = plot(x, f, 'k');
% label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'south', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Pre-call exp onset (ms)') % doesn't look right.... it's negative and too long
ylabel('Post-call insp onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

% look only at calls that occur during expirations
expInd = find(x < -250);
x = x(expInd);
y = y(expInd);

% remove outliers
[x_noOutliers, outlierIdx] = rmoutliers(x);
y_noOutliers = y(find(outlierIdx == 0));

[r, pval] = corr(x_noOutliers', y_noOutliers') % weak negative correlation, statistically significant
[p, S] = polyfit(x_noOutliers,y_noOutliers,1);
[f, delta] = polyval(p, x_noOutliers, S');

scatter(abs(x_noOutliers), y_noOutliers, 'r')
h = plot(abs(x_noOutliers), f, 'r');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'west', 'fontsize', 11, 'fontweight', 'bold')
xlabel('Pre-call exp onset (ms)') % doesn't look right.... it's negative and too long
ylabel('Post-call insp onset (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)

% look only at calls that occur during inspirations
x = (expOn.preCall1 - callWin_pre) * 1000 / fs;
y = (inspOn.postCall1 - callWin_pre) * 1000 / fs;

expInd = find(x > -250);
x = x(expInd);
y = y(expInd);

% remove outliers
[x_noOutliers, outlierIdx] = rmoutliers(x);
y_noOutliers = y(find(outlierIdx == 0));

[r, pval] = corr(x_noOutliers', y_noOutliers') % weak negative correlation, statistically significant
[p, S] = polyfit(x_noOutliers,y_noOutliers,1);
[f, delta] = polyval(p, x_noOutliers, S');

scatter(abs(x_noOutliers), y_noOutliers, 'b')
h = plot(abs(x_noOutliers), f, 'b');
label(h, strcat('r = ', num2str(r), ' p = ', num2str(pval)), 'location', 'east', 'fontsize', 11, 'fontweight', 'bold')


% legend('outliers', 'call onset during insp', 'call onset during exp?')

%% === compare pre-call exp time to time of following first exp ===
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
[~, inspDur_preCall_sortInd] = sort(inspDura.preCall1);
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

T = linspace(-preWin,postWin,length(callArray_sorted(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing normalized to max inspiration')
xlabel('Time relative to call onset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% === sort wrt pre-call expiration onset ===
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

%% === sort wrt post-call insp onset ===
[~, inspOn_postCall_sortInd] = sort(inspOn.postCall1);
callArray_sorted = callArray(:, inspOn_postCall_sortInd);
callArray_sorted_norm = (callArray_sorted ./ -min(callArray_sorted));
y = (1 : 1 : length(inspOn.postCall1))';

figure
imagesc(callArray_sorted_norm')
hold on
scatter(inspOn.preCall1(inspOn_postCall_sortInd), y, 'w', 'filled')
scatter(inspOn.postCall1(inspOn_postCall_sortInd), y, 'w', 'filled')
scatter(expOn.preCall1(inspOn_postCall_sortInd), y, 'r', 'filled')
scatter(expOn.postCall1(inspOn_postCall_sortInd), y, 'r', 'filled')

T = linspace(-preWin,postWin,length(callArray_sorted(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing normalized to max inspiration')
xlabel('Time relative to call onset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% === sort wrt call Vt ===
 % absVt vs expVt.call
[~, expVt_sortInd] = sort(expVt.call);
callArray_sorted = callArray(:, expVt_sortInd);
callArray_sorted_norm = (callArray_sorted ./ -min(callArray_sorted));
y = (1 : 1 : length(inspOn.postCall1))';

figure
imagesc(callArray_sorted_norm')
hold on
scatter(inspOn.preCall1(expVt_sortInd), y, 'w', 'filled')
scatter(inspOn.postCall1(expVt_sortInd), y, 'w', 'filled')
scatter(expOn.preCall1(expVt_sortInd), y, 'r', 'filled')
scatter(expOn.postCall1(expVt_sortInd), y, 'r', 'filled')

T = linspace(-preWin,postWin,length(callArray_sorted(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing sorted by call expiratory Vt')
xlabel('Time relative to call onset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% === sort wrt preCall Vt ===
[~, absVt_sortInd] = sort(absVt);
callArray_sorted = callArray(:, absVt_sortInd);
callArray_sorted_norm = (callArray_sorted ./ -min(callArray_sorted));
y = (1 : 1 : length(inspOn.postCall1))';

figure
imagesc(callArray_sorted_norm')
hold on
scatter(inspOn.preCall1(absVt_sortInd), y, 'w', 'filled')
scatter(inspOn.postCall1(absVt_sortInd), y, 'w', 'filled')
scatter(expOn.preCall1(absVt_sortInd), y, 'r', 'filled')
scatter(expOn.postCall1(absVt_sortInd), y, 'r', 'filled')

T = linspace(-preWin,postWin,length(callArray_sorted(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing sorted by pre-call absolute Vt')
xlabel('Time relative to call onset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% === sort wrt pre-call Vt + call Vt ===
[~, absCallVt_sortInd] = sort(absVt + expVt.call);
callArray_sorted = callArray(:, absCallVt_sortInd);
callArray_sorted_norm = (callArray_sorted ./ -min(callArray_sorted));
y = (1 : 1 : length(inspOn.postCall1))';

figure
imagesc(callArray_sorted_norm')
hold on
scatter(inspOn.preCall1(absCallVt_sortInd), y, 'w', 'filled')
scatter(inspOn.postCall1(absCallVt_sortInd), y, 'w', 'filled')
scatter(expOn.preCall1(absCallVt_sortInd), y, 'r', 'filled')
scatter(expOn.postCall1(absCallVt_sortInd), y, 'r', 'filled')

T = linspace(-preWin,postWin,length(callArray_sorted(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing sorted by pre-call absolute Vt + call expiratory Vt')
xlabel('Time relative to call onset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% === align to call expiration offset ===
callArray_postCalloffset = [];
for i = 1 : length(inspOn.postCall1)
    callArray_postCalloffset = [callArray_postCalloffset callArray(inspOn.postCall1(i) - expOn.call(i) : end - (max(inspOn.postCall1) - inspOn.postCall1(i))', i)];
end

callArray_postCalloffset_norm = (callArray_postCalloffset ./ -min(callArray_postCalloffset));
y = (1 : 1 : length(inspOn.postCall1))';

% === not sorted ===
figure
imagesc(callArray_postCalloffset_norm')
hold on

T = linspace(-(inspOn.postCall1(i) - expOn.call(i)) * 1000 / fs , (length(callArray_postCalloffset(:, 1)) - (inspOn.postCall1(i) - expOn.call(i))) * 1000 / fs,length(callArray_postCalloffset(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing aligned to call expiratory offset - not sorted')
xlabel('Time relative to call offset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% sort wrt post call exp duration offset
[~, expOn_preCall_sortInd] = sort(inspOn.postCall1);
callArray_postCalloffset_sorted = callArray_postCalloffset(:, expOn_preCall_sortInd);
callArray_postCalloffset_sorted_norm = (callArray_postCalloffset_sorted ./ -min(callArray_postCalloffset_sorted));
y = (1 : 1 : length(inspOn.postCall1))';


figure
imagesc(callArray_postCalloffset_sorted_norm')
hold on

T = linspace(-(inspOn.postCall1(i) - expOn.call(i)) * 1000 / fs , (length(callArray_postCalloffset(:, 1)) - (inspOn.postCall1(i) - expOn.call(i))) * 1000 / fs,length(callArray_postCalloffset(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing aligned to call expiratory offset - sorted by post call inspiration onset')
xlabel('Time relative to call offset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% sort wrt pre call exp onset
[~, expOn_preCall_sortInd] = sort(expOn.preCall1);
callArray_postCalloffset_sorted = callArray_postCalloffset(:, expOn_preCall_sortInd);
callArray_postCalloffset_sorted_norm = (callArray_postCalloffset_sorted ./ -min(callArray_postCalloffset_sorted));
y = (1 : 1 : length(inspOn.postCall1))';

figure
imagesc(callArray_postCalloffset_sorted_norm')
hold on

T = linspace(-(inspOn.postCall1(i) - expOn.call(i)) * 1000 / fs , (length(callArray_postCalloffset(:, 1)) - (inspOn.postCall1(i) - expOn.call(i))) * 1000 / fs,length(callArray_postCalloffset(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing aligned to call expiratory offset - sorted by preCall expiration')
xlabel('Time relative to call offset (ms)')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
colormap('jet')
colorbar

%% sort wrt pre call insp onset
[~, inspOn_preCall_sortInd] = sort(inspOn.preCall1);
callArray_postCalloffset_sorted = callArray_postCalloffset(:, inspOn_preCall_sortInd);
callArray_postCalloffset_sorted_norm = (callArray_postCalloffset_sorted ./ -min(callArray_postCalloffset_sorted));
y = (1 : 1 : length(inspOn.postCall1))';

figure
imagesc(callArray_postCalloffset_sorted_norm')
hold on

T = linspace(-(inspOn.postCall1(i) - expOn.call(i)) * 1000 / fs , (length(callArray_postCalloffset(:, 1)) - (inspOn.postCall1(i) - expOn.call(i))) * 1000 / fs,length(callArray_postCalloffset(:, 1)));
xticks = 1 : (length(T) - 1) / 6 : length(T);
xticklabels = linspace(-preWin, postWin, numel(xticks));
title('Breathing aligned to call expiratory offset - sorted by preCall inspiration')
xlabel('Time relative to call offset (ms)')
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


%% === call rhythmicity ===

% === inter-call interval distribution ===
figure
subplot(2, 1, 1)
nhist(callParams.interCallinterval.single * 1000 / fs)
title('Inter-Call Interval')
xlabel('Time (ms)')

subplot(2, 1, 2)
nhist(log(callParams.interCallinterval.single * 1000 / fs))
title('Inter-Call Interval with log')
xlabel('log_1_0 Time (ms)')
set(gca, 'tickdir', 'out', 'fontsize', 15)


figure
subplot(2, 1, 2)
nhist(log(callParams.interCallinterval.all * 1000 / fs)) % is this correct? why does it look weird?
xlabel('log_1_0  Time (ms)')
title('Call time autocorrelation with log')

subplot(2, 1, 1)
nhist(callParams.interCallinterval.all * 1000 / fs)
xlabel('Time between calls in sequence (ms)')
title('Call time autocorrelation')


% === breath autocorrelation ===



% === call autocorrelation ===

%% === call rate vs respiratory rate ===


% === look at how resp phase compares to song onset and how it's effected

for i = 1 : 5
    figure; hold on
    a(1) = subplot(2, 1, 1);
    plot(data{i}.song, 'k')
    a(2) = subplot(2, 1, 2);
    plot(data{i}.air, 'k')
    linkaxes(a, 'x')
end






























