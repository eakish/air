%% === dynamical systems analysis of air calling ===

load('callsForPhaseAnaly.mat') % load data for phase plane analysis


%% === phase plots === 
fs = 32000;
inspAmpThresh = -500; expAmpThresh = 400; durThresh = 2000;
bird = 'or82or24'; callThresh = 5000;
downSampTimes = 100;

% === iterate through data ===
for i = 1 :length([forHamish.callOnset])
    if isempty(forHamish(i).breathing)
        continue
    end
    
% % === downsample data ===
%     V = downsample(forHamish(i).breathing, downSampTimes);
% 
% % === calculate the derivative ===
%     dv =  diff(V);
% 
% % % === plot V vs dV/dt ===
% %     figure(1)
% %     plot(V(1 : end - 1), dv); hold on
% 
%     figure
%     subplot(2, 1, 1)
%     plot(V(1 : end - 1), dv)
%     xlabel('V')
%     ylabel('dV/dt')
%     set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15) 
%     
%     subplot(2, 1, 2)
%     plot(V(1 : end - 1))
%     xlabel('Time')
%     ylabel('V')
%     set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)


    
    % find expiratory crossings
    [breathStartInd, breathEndInd, breathRate,  inspVol, expVol, inspEnd, inspPeak, expPeak, inspDur, expDur] = ek_segmentBreaths_edit(forHamish(i).breathing, bird, durThresh, inspAmpThresh, expAmpThresh);
    
    % find expiratory period
    breathPeriod = diff(inspEnd) * 1000 / fs;
    
    % identify which one has a call - based on crossing threshold
    callID = find(expPeak > callThresh);
    
    if callID > length(breathPeriod)
        continue
    end
    
    isCall = zeros(length(breathPeriod) - 1, 3);
    isCall(callID - 1, :) = [1 0 0]; % labels n before call and n + 1 containing the call
    
%     figure; plot(forHamish(i).breathing)
%     hold on; scatter(inspEnd, zeros(length(inspEnd), 1))
%     
%     % ======== plot poincare ========
%     figure(2)
%     scatter(breathPeriod(1 : end - 1), breathPeriod(2 : end), 50, isCall)
%     hold on
    
    
    breaths = zeros(max(breathPeriod) * fs / 1000, length(breathStartInd));
    % == create breath array to look at period ===
    for j = 1 : length(breathStartInd) - 2
        breaths(1 : length (forHamish(i).breathing(inspEnd(j) : inspEnd(j + 2))), j) = forHamish(i).breathing(inspEnd(j) : inspEnd(j + 2Z));
    end
    
    
    preWin = 0; postWin = 0;
    ek_mouseOverBreaths('period n', breathPeriod(1 : end - 1), 'period n + 1', breathPeriod(2 : end), breaths, preWin, postWin, isCall)
    unity = 1 : 700;
    plot(unity, unity, 'k')


% % ====== for looking at individual call phase plots =====
% % === downsample data ===
%     V = downsample(forHamish(i).breathing(breathStartInd(callID - 2) : breathStartInd(callID + 1)), downSampTimes);
%     V_call = downsample(forHamish(i).breathing(breathStartInd(callID) : breathEndInd(callID)), downSampTimes);
% 
% %     V_call = V(
% 
% % === calculate the derivative ===
%     dv =  diff(V);
%     dv_call = diff(V_call);
% 
% % % === plot V vs dV/dt ===
% %     figure(1)
% %     plot(V(1 : end - 1), dv); hold on
% 
%     figure
%     subplot(2, 1, 1)
%     plot(V(1 : end - 1), dv, 'k'); hold on
%     plot(V_call(1 : end - 1), dv_call, 'r')
%     xlabel('V')
%     ylabel('dV/dt')
%     title(strcat('downsampled x', num2str(downSampTimes)));
%     set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15) 
%     
%     x = linspace(0, length(V(1 : end - 1)) * 1000 / fs * downSampTimes, length(V(1 : end - 1)));
% %     x_call = linspace(1, length(V_call(1 : end - 1)) * 1000 / fs * downSampTimes, length(V_call(1 : end - 1)));
%     subplot(2, 1, 2)
%     plot(x, V(1 : end - 1), 'k'); hold on
% %     plot(V_call(1 : end - 1), 'k')
%     xlabel('Time (ms)')
%     ylabel('V')
%     set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)



end

% figure(1)
% xlabel('V')
% ylabel('dV/dt')
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)
%

% figure(2)
% unity = 1 : 700;
% plot(unity, unity, 'k')
% xlabel('Expiratory Period n')
% ylabel('Expiratory Period n + 1')
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'tickdir', 'out', 'fontsize', 15)

% look at when the respiratory trace diverges

















