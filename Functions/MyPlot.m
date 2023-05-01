%% My Parameters of Figure Plotting

%% Function
function MyPlot(fplot,resize,multi)
    if multi == 0
        % Create a Window at a Decided Position
        figure('Color','white','Position',[150,70,(150+950*resize),(70+830*resize)]);
        % Create an Frame for a Figure
        axes('OuterPosition',[0.005,0.01,.99,.97]);
    elseif multi == 1
        clear sub;
        % Open the Overwrite Mode
        hold on
    end
    % Create a Specific Plot
    plot(fplot.x,fplot.y,fplot.colorstyle,'linewidth',2);
    % Set log Scale
    if isfield(fplot,'yscale') == 1
        set(gca,'YScale', fplot.yscale);
    end
    % Box Setting
    box on
    % Box Linewidth
    set(gca,'linewidth',2);
    % Grid Setting
    grid off
    % Axis Range
    axis(fplot.range);
    % Font Size
    set(gca,'fontsize',30*resize,'fontname','Times New Roman');
    % Label of x Axis
    xlabel(fplot.xlabel,'interpreter','latex');
    % Label of y Axis
    ylabel(fplot.ylabel,'interpreter','latex');
    % Figure Title
    if isfield(fplot,'title') == 1
        title(fplot.title);
    end
end
