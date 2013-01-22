function plotribbon(startPosition,stepDistance,endPosition,superElPlot,myLegend,myxticklabel,outpath,newfilename,saveString)       

    figure1 = figure
    axes1 = axes('Parent',figure1,'Projection','perspective','FontSize',12);
    x = startPosition:stepDistance:endPosition;
    %hold all;
    %x = (1:nSteps+1)*;
    %ylabel('t (s)','FontSize',12)
    ribbon(x,transpose(superElPlot))
    %ribbon(x,transpose(superElPlot(21,:)))
    view([-66.5 26]);
    ylabel('Axial ion position (mm)','FontSize',12)
    xlabel('Electrode Number','FontSize',12)
    zlabel('Voltage (V)','FontSize',12)
    xlim([1 21])
    zlim([min(min(superElPlot)) max(max(superElPlot))])
    ylim([min(x) max(x)])
    set(gca,'xtick',[1 5 10 15 18 21])
    set(gca,'xticklabel',myxticklabel)
    legend1 = legend(myLegend);
    set(legend1,'FontSize',9)
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperUnits','inches','PaperPosition',[-0.35 -0.2 8 6])
    set(gcf, 'PaperSize', [7 5.4]);
    print(figure1,'-dpdf',[sprintf('%s%s/pdf/%s.pdf',outpath,newfilename,saveString)]) 
    saveas(figure1,[sprintf('%s%s/fig/%s.fig',outpath,newfilename,saveString)])

    min(min(superElPlot))
    max(max(superElPlot))
    
    %hold off;