function shadedTimeSeries(X, Ydata, indicator, Xlabel, Ylabels, colour, Yspacing)
% SHADEDTIMESERIES(X, Ydata, indicator, Xlabel, Ylabels, colour, Yspacing)
%
% Plot time series one above the other with coloured strips highlighting
% interesting features.
%
% All data is in columns.
%
% X is the time vector.
% Each column of Ydata will be plotted in a separate subplot, one beneath
% the other.
% Indicator is a vector of ones and zeroes where the ones indicate the
% areas of interest that will be shaded.
% Xlabel is a string for labeling the bottom horizontal axis.
% Ylabels is a cell array of strings used for labeling the Y axes.
% Colour is an RGB colour array that will be used for shading.
% Yspacing is the percentage of vertical padding to leave above and below
% each plot so that they don't touch the top and bottom of their respective
% plots.
%
% Example: 
%           time = [1:50]';
%           acceleration = rand(50,1);
%           velocity = cumtrapz(acceleration, time);
%           position = cumtrapz(velocity, time);
%           indicator = [zeros(1,10) ones(1,10) zeros(1,25) ones(1,5)]';
%           shadedTimeSeries(time, [acceleration velocity position], indicator, 'Time', {'Acceleration' 'Velocity' 'Position'}, [1 .7 1], 10);
%
% Carl Fischer. July 2008.
% http://eis.comp.lancs.ac.uk/~carl/blog/

if nargin < 3
    error('Please provide at least X, Ydata and indicator parameters.');
end
if nargin < 4
    Xlabel = '';
end
if nargin < 5
    Ylabels = cell(1,size(Ydata,2));
end
if nargin < 6
    colour = [1 .7 .7];
end
if nargin < 7
    Yspacing = 5;
end



start_marks = find(diff(indicator) > 0);
end_marks = find(diff(indicator) < 0);
if start_marks(1) > end_marks(1) % plot is shaded from start
    start_marks = [1; start_marks];
end        
if start_marks(end) > end_marks(end) % plot is shaded until end
    end_marks = [end_marks; size(X,1)];
end

figure;
hold on
for i = 1:size(Ydata,2)
    plot_partial(X, Ydata(:,i), i, Ylabels{i}, size(Ydata,2));
    if i ~= size(Ydata,2)
        set(gca,'xticklabel',[]) % remove x label for all except bottom plot
    end
end
xlabel(Xlabel);



    function plot_partial(x, ydata, number, name, total)
        subplot(total,1,number);
        hold on
        ylabel(name);
        xlim([x(1) x(end)]);% force plot to go from edge to edge
        padding = Yspacing/100*(max(ydata)-min(ydata));% spacing at top and bottom of plotted line
        ylim([min(ydata)-padding max(ydata)+padding]); % default leaves too much space around graph
        
        for j = 1:min(size(start_marks,1), size(end_marks,1)) % create all the shaded areas
            % This patch line was ripped from ShadePlotForEmphasis.m by
            % Michael Robbins on Matlab Central.
            patch([repmat( x(start_marks(j)),1,2) repmat( x(end_marks(j)),1,2)], ...
                [get(gca,'YLim') fliplr(get(gca,'YLim'))], ...
                [0 0 0 0],colour,'EdgeColor', 'none');
        end
        plot(x, ydata, 'LineWidth', 3); % plot data line on top of colour patch
        set(gca, 'layer', 'top'); % put ticks back on top of colour patch
    end
end
