
% Created by EK in April 2020
% do PCA on call breathing array to see if there's clustering of call types


% rd8Calls = load('rd8_110719_airCallsArray.mat');
orCalls = load('or84or24_121819_airCallsArray.mat');

%%
% breathArray_calls = [rd8Calls.breathArray_calls orCalls.breathArray_calls];
% whichBird = zeros(length(breathArray_calls(1, :)), 3);
% whichBird(1 : length(rd8Calls.breathArray_calls(1, :)), 1) = ones(length(rd8Calls.breathArray_calls(1, :)), 1);
%% === cut out section around call ===
callMat = breathArray_calls(2000 * fs / 1000 + 100 : 67500, :);

% callMat = smoothdata(callMat, 'SmoothingFactor', 0.05);

figure; plot(callMat)

%% === do low-dimensional mapping on breathing traces ===
coeff = pca(callMat); % PCA
%
mapping_tsne = tsne(callMat'); % TSNE

%% === plot calls in PCA space to see if there's clustering of calls ===
figure; subplot(2, 1, 1); scatter(coeff(:, 1), coeff(:, 2))
title('call PCA')


subplot(2, 1, 2); scatter(mapping_tsne(:, 1), mapping_tsne(:, 2))
title('call TSNE')

%% === look at what the waveforms of these calls looks like for different clusters ===
ek_mouseOverBreaths('tsne 1', mapping_tsne(:, 1), 'tsne 2', mapping_tsne(:, 2), callMat, preWin, postWin) %, whichBird)

%% for PCA

ek_mouseOverBreaths('PCA 1', coeff(:, 1), 'PCA 2', coeff(:, 2), callMat, preWin, postWin) %, whichBird)


%% === map this onto the actual vocalizations ===







%% === do the relationships between respiratory phase and calls hold for all call types? ===








































