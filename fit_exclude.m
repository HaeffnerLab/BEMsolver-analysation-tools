function out = fit_exclude(startInt,endInt,stepInt,superElData) 

%startInt = 0;
%endInt = 0.1;
%stepInt = 0.01;

x = startInt:stepInt:endInt;


figure; 
hold all;
%axes1 = axes('Parent',figure,'DataAspectRatio',[5 2.5 1]);
%box(axes1,'on');
ylabel({'Voltage [V]'});
xlabel({'z-axis [mm]'});

grid on;
box on;


for iii=1:21

    x = x(:);
    y1 = superElData(iii,:)';

    ex_ = false(length(x),1);

    ex2_=[];
    
    
    % check distances between datapoints: if there is a discontinuity,
    % exclude the points between jump up and down
    
    for ii=1:1:size(y1,1)
       if ~(ii==1)
           if (abs(y1(ii)-y1(ii-1))>0.3*abs(y1(ii)))            
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

    plot(x,cf_(x))
    
    superCf(:,iii) = cf_(startInt:stepInt:endInt);

end

legend1 = legend('W1 ','W2 ','W3 ','W4 ','W5 ','W6 ','W7 ','W8 ','W9 ','W10','N1 ','N2 ','N3 ','N4 ','N5 ','N6 ','N7 ','N8 ','N9 ','N10','CNT','RF ');
set(legend1,'Location','EastOutside');
hold off;

out=superCf;