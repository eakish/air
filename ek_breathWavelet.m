
approx_mat = []; % for approximation coefficients
cd1_mat = [];
for i = 1 : length(callMat(1, :))
    [c,l] = wavedec(callMat(:, i),4,'db2');
    approx = appcoef(c,l,'db2');
    [cd1,cd2,cd3] = detcoef(c,l,[1 2 3]);
    
    % remove weird recordings -- ones that are noisy or something. should do this earlier.
    if min(approx) < 1000
        continue
    end
    
    approx_mat = [approx_mat approx];
    cd1_mat = [cd1_mat cd1];
% 
%     % Plot the coefficients.
%     figure
%     subplot(4,1,1)
%     plot(approx)
%     title('Approximation Coefficients')
%     subplot(4,1,2)
%     plot(cd3)
%     title('Level 3 Detail Coefficients')
%     subplot(4,1,3)
%     plot(cd2)
%     title('Level 2 Detail Coefficients')
%     subplot(4,1,4)
%     plot(cd1)
%     title('Level 1 Detail Coefficients')
    
end
%
figure; plot(approx_mat)




%% pca on approximation coefficients
coeff_approx = pca(approx_mat(150 : end, :));

k_approx = kmeans(coeff_approx, 2); % this is definitely not working
c = [k_approx k_approx k_approx];
c = c - 1;
x = find(c(:, 1) == 1);
c(x, 2) = 0;
c(x, 3) = 0;
% coeff_cd1 = pca(cd1_mat);


% plot pca coefficients
%%
figure; scatter(coeff_approx(:, 1), coeff_approx(:, 2), 50, c)

% figure; scatter(coeff_cd1(:, 1), coeff_cd1(:, 2))

%% implement zhang and ghazanfar analysis?




%% k-means clustering?







