% gui that displays respiratory waveform after clicks

function ek_mouseOverBreaths(xLabel, xArray, yLabel, yArray, breaths, preWin, postWin, category) %, breathDur, inspPeak, expPeak)


    % == click on a data point and display the original curve ==
    figure
    subplot(211)
%     h = scatter(xArray, yArray, 50, 'k');
    h = scatter(xArray, yArray, 50, category);
    xlabel(xLabel)
    ylabel(yLabel)
    unity = min(min(xArray, yArray)) : max(max(xArray, yArray));
    hold on
%     plot(unity, 'r')
    set(gca, 'tickdir', 'out', 'fontsize', 15)

    set(h,'buttondownfcn', {@hitme, breaths, xArray, yArray, preWin, postWin});

    %% GUI
    function hitme(gcbo, ~, breaths, xArray, yArray, preWin, postWin)

    [x,y] = ginput(1);

    cpFigure = get(gcf, 'CurrentPoint');
    cpAxis = get(gca, 'CurrentPoint');

    % Retrieve the x and y data from the plot
    xdata = get(gcbo, 'xdata');
    ydata = get(gcbo, 'ydata');

    % Scan the actual plotted points, figuring out which one comes closest
    distances = sqrt((x-xdata).^2+(y-ydata).^2);
    [minValue minIndex] = min(distances);

    subplot (211)
    hold on;
    scatter(xArray(minIndex), yArray(minIndex), 'filled')

    subplot (212)
%     x = linspace(preWin, postWin, length(breaths(:, minIndex)));
    plot(breaths(:, minIndex)); hold on
    plot(zeros(length(x), 1), 'k')
    
    xlabel('time from call onset')
    set(gca, 'tickdir', 'out', 'fontsize', 15)
    ylim([-5000 15000])
    
    end

end