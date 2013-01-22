function imagesc_plot(superData,superEl,outpath,newfilename,nSteps,stepDistU,startPosU,endPosU)

    close all

    RFampl = 150; 
    Freq = 40e6;
    e=1.60217646e-19;                                                                       % elementary charge unit in SI
    mp=1.67262158e-27;                    
    Omega = 2*pi*Freq;

    h=figure;

    
    l = size(superEl,2);
    
    vec = round(1/2*(cos(pi:pi/15:2*pi)+1)*l);
    vec = vec(2:size(vec,2));
    movePosVec = vec
    
    %if -630 to 630, if truncated comment out!!
    %movePosVec=[1,3,14,23,33,48,63,78,93,102,113,124,126];
    
    for i=(1:size(movePosVec,2))
        movePos=movePosVec(i);

        %        actPos = startPosU + stepDistU * (movePos+14);

        actPos = startPosU + stepDistU * (movePos+(127-l)/2);
        
        clear newData;
        newData.Z = [];
        Unew = [];

        for pos = 0:nSteps

            data = superData(pos+1);

                Vl=CalcVDC(data,superEl(:,movePos),0,0,0,40,2,0,0,0,0);
                X = data.X; Y =data.Y; Z =data.Z;
                [y x z] = meshgrid(Y,X,Z);
                grid = data.grid;
                Vrf = RFampl*data.EL_RF;
                [Ex,Ey,Ez] = gradient(Vrf,1e-3*grid(4),1e-3*grid(5),1e-3*grid(6));
                Esq1 = Ex.^2 + Ey.^2 + Ez.^2;
                PseudoPhi = e^2*Esq1/(4*40*mp*Omega^2);
                U = PseudoPhi+e*Vl;

                aa = U(:,:,1:10)/e;
                aaa = permute(aa,[3 1 2]);

                Unew = vertcat(Unew,aaa);

        end
        


        subplot('position',[0.2 0.91-(0.06)*(i-1) 0.6 0.04])

        hold all;
        box on;
        set(gca,'xtick',[])
        set(gca,'FontSize',12)

        if (i==round(size(movePosVec,2)/2)) 
            ylabel('x (mm)','FontSize',12);
        end

        x = ((0:nSteps*stepDistU)-nSteps*stepDistU/2)*0.001;
        y = ((0:2:40)-20)*0.001;

        ylim([-0.020 0.020]);
        xlim([startPosU/1000 endPosU/1000]);

        i 
        %,[-1.5 -1]
        imagesc(x,y,Unew(:,:,10)',[-1.5 -1])
        
        plot(actPos/1000,0,'o','MarkerFaceColor','w','MarkerSize',5,'MarkerEdgeColor','none')
        

        hold off;
    end

    set(gca,'xtickMode', 'auto','FontSize',12)
    xlabel('z (mm)','FontSize',12);
    %colorbar([0.83 0.85-(0.18-0.035)*(i-1) 0.03 0.85-(0.18-0.035)*(1)-(0.85-(0.18-0.035)*(i))+0.1],'FontSize',12);

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperUnits','inches','PaperPosition',[0.07*1.1 -0.15*1.1 3.8*1.1 8.5*1.1])
    set(gcf, 'PaperSize', [3.4*1.1 8.2*1.1]);

    print(h,'-dpdf',[sprintf('%s%s/',outpath,newfilename) 'imagesc_potential_traj.pdf'])
    saveas(h,[sprintf('%s%s/fig/',outpath,newfilename) 'imagesc_potential_traj.fig'])





