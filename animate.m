function animate(x, y, rate, limits)

    p = plot(x, squeeze(y(:,1,:)));
    xlim([min(x), max(x)]);
    if nargin == 4
        ylim(limits);
    else
        ylim([min(y, [], "all"), max(y, [], "all")]);
    end

    m = size(y, 2);
    n = size(y, 3);

    for i = 1:rate:m
        try
            for j = 1:n
                p(j).YData = y(:,i,j);
            end
        catch e
            if e.identifier == "MATLAB:class:InvalidHandle"
                warning("Not able to update plot, " ...
                    + "assumed that user ended early.");
                break
            else
                rethrow(e);
            end
        end
        drawnow;
    end

end
