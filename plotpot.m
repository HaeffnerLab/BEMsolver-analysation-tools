function plotpot(A,I,J,K,grid,f,tit,ylab,outpath,newfilename)
% plotpot(A,I,J,K,grid,f,tit,ylab)
% Make 2d mesh plots and 1d plots of A around the point I,J,K
% f==0: no plots, 1: 2d plots, 2: 1d plots, 3: both 2d and 1d plots 
% remember: I->X, J->Y, K->Z
% tit is the title of the plots produced
% ylab is the label on the y axis of the plot produced

if ~f, return; end;

if ( (f==1)||(f==3) ),
    meshslice(A,1,grid);
    meshslice(A,2,grid);
    meshslice(A,3,grid);
end

if ( (f==2)||(f==3) ),
    % plot I
    PrI = A(:,J,K); 
    t = grid(1):grid(4):grid(1)+(size(A,1)-1)*grid(4);
    close all; 
    h = figure;
    axes1=axes('Parent',h,'FontSize',12);
    box(axes1,'on');
    hold(axes1,'all');
    plot(t,PrI); plot(t(I),PrI(I),'r*'); 
    xlabel('x (mm)','FontSize',12); ylabel(ylab,'Interpreter','tex');
   % ylim([min(PrI) max(PrI)])
    xlim([min(t) max(t)])
    %print(h,'-depsc',[sprintf('%s%s/',outpath,newfilename) 'pot_radial.eps'])
    %saveas(h,[sprintf('%s%s/',outpath,newfilename) 'pot_radial.fig'])
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0.1 3.2 2.4])
    set(gcf, 'PaperSize', [3.1 2.4]);
    %print(h,'-dpdf',[outpath 'eps/' sprintf('%i',newfilename) sprintf('_%s_radial.pdf',tit)])
    %saveas(h,[outpath 'fig/' sprintf('%i',newfilename) sprintf('_%s_radial.fig',tit)])
    hold off; 
    pause;
    
    % plot J
    PrJ = A(I,:,K); 
    t = grid(2):grid(5):grid(2)+(size(A,2)-1)*grid(5);
    close all;     h = figure;
    axes1=axes('Parent',h,'FontSize',12);  
    box(axes1,'on');
    hold(axes1,'all');
    plot(t,PrJ); plot(t(J),PrJ(J),'r*');
    xlabel('y (mm)','FontSize',12); ylabel(ylab,'Interpreter','tex');
    %ylim([min(PrJ) max(PrJ)])
    xlim([min(t) max(t)])
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0.1 3.2 2.4])
    set(gcf, 'PaperSize', [3.1 2.4]);
    %print(h,'-dpdf',[outpath 'eps/' sprintf('%i',newfilename) sprintf('_%s_height.pdf',tit)])
    %saveas(h,[outpath 'fig/' sprintf('%i',newfilename) sprintf('_%s_height.fig',tit)])
    hold off; 
    pause;
    
    % plot K
    aaa = permute(A,[3 1 2]);
    PrK = aaa(:,I,J);
    t = grid(3):grid(6):grid(3)+(size(A,3)-1)*grid(6);
    close all;     
    h = figure;
    axes1=axes('Parent',h,'FontSize',12);
    box(axes1,'on');
    hold(axes1,'all');
    plot(t,PrK); plot(t(K),PrK(K),'r*');
    xlabel('z (mm)','FontSize',12); ylabel(ylab,'Interpreter','tex'); 
    %ylim([min(PrK) max(PrK)])
    xlim([min(t) max(t)])
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0.1 3.2 2.4])
    set(gcf, 'PaperSize', [3.1 2.4]);
    %print(h,'-dpdf',[outpath 'eps/' sprintf('%i',newfilename) sprintf('_%s_axial.pdf',tit)])
    %saveas(h,[outpath 'fig/' sprintf('%i',newfilename) sprintf('_%s_axial.fig',tit)])
    hold off; 
    pause;
    
    %pause;
end