function imagesc_plot(superData,totv,nSteps,stepDistU,startPosU,endPosU)

    RFampl = 150; 
    Freq = 40e6;
    e=1.60217646e-19;                                                                       % elementary charge unit in SI
    mp=1.67262158e-27;                    
    Omega = 2*pi*Freq;

    h=figure;

   
    
    
              
        clear newData;
        newData.Z = [];
        Unew = []; 

        for pos = 0:nSteps

            data = superData(pos+1);

                Vl=CalcVDC(data,totv,0,0,0,40,2,0,0,0,0);
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
        
        
        %hold all;
        box on;
        set(gca,'xtick',[])
        set(gca,'FontSize',12)


        x = ((0:nSteps*stepDistU)-nSteps*stepDistU/2)*0.001;
        y = ((0:2:40)-20)*0.001;

        
        %ylim([-2.0 -0.8]);
        
        xlim([startPosU/1000 endPosU/1000])

        plotUnew = Unew(1:1261,10,10);
        
        
        
        plot(x,plotUnew)
        pause;
        imagesc(x,y,Unew(:,:,10)')
                        

        
    set(gca,'xtickMode', 'auto','FontSize',12)
    xlabel('z (mm)','FontSize',12);

  




