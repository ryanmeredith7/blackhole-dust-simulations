% This function shows an animation of one or more data sets. Can be used to
% animate more than one data set on top of each other by calling it as
% animate(x, cat(3, y1, y2, ..., yn), rate).
function animate(x, y, rate, limits)

    arguments

        % Input data for the x axis
        x (:,1) {mustBeReal}

        % Input data for the Y axis
        y (:,:,:) {mustBeReal}

        % Rate at which the animation plays.
        rate (1,1) {mustBePositive, mustBeInteger}

        % Optional input for the limits of the y axis
        limits (1,2) {mustBeReal} = [min(y, [], "all"), max(y, [], "all")]

    end

    % Draws the initial plot
    p = plot(x, squeeze(y(:,1,:)));
    xlim([min(x), max(x)]);
    ylim(limits);

    % Gets the length of the animation and the number of data sets to plot
    m = size(y, 2);
    n = size(y, 3);

    % Does the animation: updates the plots if the user didn't exit the plot,
    % and redraws the plot
    for i = 1:rate:m
        try
            for j = 1:n
                p(j).YData = y(:,i,j);
            end
        catch e
            if e.identifier == "MATLAB:class:InvalidHandle"
                warning("Animate:FailedToUpdatePlot", ...
                    "Not able to update plot, assumed that user ended early.");
                break
            else
                rethrow(e);
            end
        end
        drawnow;
    end

end
