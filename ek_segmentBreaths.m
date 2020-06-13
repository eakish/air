function [breathStartInd, breathEndInd, breathRate,  inspVt, expVt, inspEnd, inspPeak, expPeak, inspDur, expDur] = ek_segmentBreaths(breathArray, bird, durThresh, inspAmpThresh, expAmpThresh)
% === Function for segmenting voltage traces into individual breaths ---
% somethings wrong with the part that splits it into inspirations and
% expirations
tic;

% if it crosses above/below zero for more than x ms, then that's start/stop of breath
% save indices of breath starts and stops (can also use this to find actual time)

% bird = 'cage106_10';
% experiments = fieldnames(data.(bird));
% durThresh = 100;
% inspAmpThresh = -500;
% expAmpThresh = 0;

legend_conds = {};
% bird = 'gr98bu75';
% condition = 'bird';
c = 'k';
for iter = 1 : 1 %length(experiments)

    crossZero =[];
    threshCnt = 0; % count up time from now until thresh - account for noise/minibreaths
    isBreath = 0;
    breathStart = 1; % idx of start (col1) and end (col2) of each breath (row)
    breathEndInd = [];
    crossZero = 1;
    exp = 0;

    for i = 2 : length(breathArray) % just find crossings first
        prev = sign(breathArray(i - 1));
        curr = sign(breathArray(i));
        if prev == -1 && curr == +1
            exp = 1;
        end
        if prev == +1 && curr == -1     
                if exp % checks to make sure that an expiration happened between inspirations
                    breathStart = [breathStart i];
                    breathEndInd = [breathEndInd (i - 1)];
                    exp = 0;
                    if breathEndInd(end) - breathStart(end - 1) > durThresh
                        if max(breathArray(breathStart(end - 1) : breathEndInd(end))) > expAmpThresh
                            if min(breathArray(breathStart(end - 1) : breathEndInd(end))) > inspAmpThresh || length([find(breathArray(breathStart(end - 1) : breathEndInd(end)) < 0)]) < 500
                                breathStart = [breathStart(1 : end - 2) i];
%                                 breathEndInd = [breathEndInd(1 : end - 1) i - 1];
                                breathEndInd = [breathEndInd(1 : end - 2) i];
                                exp = 1;
                            end
                        else
                            exp = 0;
                        end
                    else
                        breathStart = [breathStart(1 : end - 2) i];
%                         breathEndInd = [breathEndInd(1 : end - 1) i - 1];
                        breathEndInd = [breathEndInd(1 : end - 2) i];
                        exp = 1;
                    end
                end
        end
    end
    breathStartInd  = breathStart(1 : end - 1); % so we cut off at the end of a breath
    breathEndInd = breathEndInd(2 : end);

%% ===== segment breaths based on indices =======
    inspPeak  = [];
    expPeak  = [];
    breathDur  = [];
    inspEnd = [];
    
    % === for calculating breath rate and tidal volume ===
    binsize = 60000; % count number of breaths per minute - currently just binning
    breathRate  = []; % breathing rate
    breathVt  = []; % tidal volume for each breath
    inspVt  = []; % inspiratory tidal volume for each breath
    expVt  = []; % expiratory tidal volume for each breath
    inspDur  = [];
    expDur  = [];
    
    breathStartInd2 = [];
    for i = 1 : length(breathStart) - 1
        breath = breathArray(breathStart(i): breathStart(i + 1) - 1);
     
        % ==== SEGMENT BREATHS INTO INSPIRATIONS AND EXPIRATIONS HERE ====
        insp = 0;
        for j = 1 : length(breath)
            if min(breath(1 : j)) < inspAmpThresh && insp == 0
                insp = 1;
            elseif insp == 1 && breath(j) >= 0
                inspEnd = [inspEnd breathStart(i) + (j - 1)];
                inspVt  = [inspVt  sum(breath(1 : j - 1))];
                expVt  = [expVt  sum(breath(j : end))];
                
                inspDur  = [inspDur  (j - 1)];
                
                break
            end
        end        
        
        if insp == 0
            continue
        end
        
        breathStartInd2 = [breathStartInd2 breathStart(i)];
        
        % split into insp and exp
        inspPeak  = [inspPeak  min(breath)];
        expPeak  = [expPeak  max(breath)];
        breathDur  = [breathDur  breathStart(i + 1) - breathStart(i)];
%         inspDur  = [inspDur  inspEnd - breathStart(i)];
        
        if isempty(inspDur)
%             figure; plot(breathArray)
%             hold on
%             scatter(breathStartInd, zeros(1, length(breathStartInd )), 'r', 'filled')
%             scatter(inspEnd, zeros(1, length(inspEnd)), 'k', 'filled')
            expDur  = [];
            break
        end
        
        expDur  = [expDur  breathDur(end) - inspDur(end)];
%         breath = breathArray(breathStart(i): breathStart(i + 1) - 1);

        
        % for calculating breath rate and tidal volume
        breathVt  = [breathVt  sum(breathArray(breathStart(i): breathStart(i + 1) - 1))];  
    end
    
    % === plot breaths and segmentation for sanity check ===
        % === plot breath array with segmented breath starts ===

    
    %% === Breathing rate across time ===
    for i = (binsize / 2) + 1 : binsize : breathStart(end) - (binsize / 2)
        breathRate  = [breathRate  length(intersect(find(breathStart >= i - (binsize / 2)), find(breathStart < i + (binsize / 2))))];
    end
%     figure(100)
%     plot(breathRate , 'LineWidth', 1.5, 'Color', c)
%     title([bird, ': ', 'Breathing rate'])
%     ylabel('Breaths Per Minute')
%     xlabel('Time (ms)')
%     set(gca,'fontsize',20,'TickDir','out')
%     hold on
%     
%     %% === Tidal volume across time ===
%     figure(200)
%     scatter(breathStartInd , breathVt , 20, c)
%     title([bird, ': ', 'Tidal Volume'])
%     ylabel('V_t (voltage)')
%     xlabel('Time (ms)')
%     set(gca,'fontsize',20,'TickDir','out')
%     hold on
    
    
    clear breathStart
end

breathEndInd = breathStartInd2(2 : end);
breathStartInd  = breathStartInd2; % so we cut off at the end of a breath

% 
% figure
% plot(breathArray)
% hold on
% scatter(breathStartInd, zeros(1, length(breathStartInd )), 'r', 'filled')
% scatter(inspEnd, zeros(1, length(inspEnd)), 'k', 'filled') % checks to make sure insp end works
































