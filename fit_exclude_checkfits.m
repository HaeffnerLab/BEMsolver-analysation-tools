startInt = 0;
endInt = 0.1;
stepInt = 0.02;

x = startInt:stepInt:endInt;


% figure; 
% hold all;
% %axes1 = axes('Parent',figure,'DataAspectRatio',[5 2.5 1]);
% %box(axes1,'on');
% ylabel({'Voltage [V]'});
% xlabel({'z-axis [mm]'});
% 
% grid on;
% box on;


for iii=1:21

    % Set up figure to receive data sets and fits
    f_ = clf;
    figure(f_);
    set(f_,'Units','Pixels','Position',[489 254 672 475]);
    % Line handles and text for the legend.
    legh_ = [];
    legt_ = {};
    % Limits of the x-axis.
    xlim_ = [Inf -Inf];
    % Axes for the plot.
    ax_ = axes;
    set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
    set(ax_,'Box','on');
    axes(ax_);
    hold on;

    % --- Plot data that was originally in data set "y1 vs. x"
        x = x(:);
        y1 = superEl(iii,:)';
    h_ = line(x,y1,'Parent',ax_,'Color',[0.333333 0 0.666667],...
        'LineStyle','none', 'LineWidth',1,...
        'Marker','.', 'MarkerSize',12);
    xlim_(1) = min(xlim_(1),min(x));
    xlim_(2) = max(xlim_(2),max(x));
    legh_(end+1) = h_;
    legt_{end+1} = 'y1 vs. x';

    % Nudge axis limits beyond data limits
    if all(isfinite(xlim_))
        xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
        set(ax_,'XLim',xlim_)
    else
        set(ax_, 'XLim',[-0.011100000000000000491, 0.1011000000000000093]);
    end

    % --- Create fit "fit 1"

    % Apply exclusion rule "tedt"
    %ex_ = (y1 >= -1);


    ex_ = false(length(x),1);

    ex2_=[];
    
    for ii=1:1:size(y1,1)
       if ~(ii==1|ii==2)
           if (abs(y1(ii)-y1(ii-1))>0.3*abs(y1(ii)))
           %if (abs(y1(ii)-y1(ii-1))>1.5*(abs(y1(ii-1))-abs(y1(ii-2))))
                 ex2_(:,ii) = ii;
           end
       end
    end

   
    for ii=min(nonzeros(ex2_)):max(nonzeros(ex2_))
       ex_(ii) = 1
    end
        

    ok_ = isfinite(x) & isfinite(y1);
    if ~all( ok_ )
        warning( 'GenerateMFile:IgnoringNansAndInfs',...
            'Ignoring NaNs and Infs in data.' );
    end
    ft_ = fittype('spline');

    % Fit this model using new data
    if sum(~ex_(ok_))<2
        % Too many points excluded.
        error( 'GenerateMFile:NotEnoughDataAfterExclusionRule',...
            'Not enough data left to fit ''%s'' after applying exclusion rule ''%s''.',...
            'fit 1', 'tedt' );
    else
        cf_ = fit(x(ok_),y1(ok_),ft_,'Exclude',ex_(ok_));
    end
    % Alternatively uncomment the following lines to use coefficients from the
    % original fit. You can use this choice to plot the original fit against new
    % data.
    %    cv_ = { 186285197576.66629028, -78549847024.221069336, 13106315211.429351807, -1075478523.7714936733, 42734435.855544321239, -613382.00978270021733, -2419.169468339197465, 134.62142721575747828, -4.4711146804558943302, -1.7650029241927982504};
    %    cf_ = cfit(ft_,cv_{:});

    % Plot this fit
    
    h_ = plot(cf_,'fit',0.95);
    set(h_(1),'Color',[1 0 0],...
        'LineStyle','-', 'LineWidth',2,...
        'Marker','none', 'MarkerSize',6);
    % Turn off legend created by plot method.
    legend off;
    % Store line handle and fit name for legend.
    legh_(end+1) = h_(1);
    legt_{end+1} = 'fit 1';

    % --- Finished fitting and plotting data. Clean up.
    hold off;
    % Display legend
    leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
    h_ = legend(ax_,legh_,legt_,leginfo_{:});
    set(h_,'Interpreter','none');
    % Remove labels from x- and y-axes.
    xlabel(ax_,'');
    ylabel(ax_,'');
    
    

    %plot(x,cf_(x))
        pause;
end

%hold off;