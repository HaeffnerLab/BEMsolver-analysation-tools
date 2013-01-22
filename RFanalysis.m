
Omega = 2*pi*30*1e6; %rf frequency in Hz
Vamp = 100; % Ampltutde Voltage Volts, the amplitude to time average conversion for the purpose of finding the pseudopotential is accounted for in the formula below
e=1.602176487e-19;                                                                       % elementary charge unit in SI
mp=1.672621637e-27;                                                            % proton mass in SI units
close all;

% L=11;
% 
% regenspacing=0.0008;
% 
% Xregen = -0.01:regenspacing:0.01;
% Yregen = 0.18:regenspacing:0.2;
% Zregen = -0.01:regenspacing:0.01;
% 
% pos = round(size(Yregen,2)/2)
% 
% 
% 
% 
% load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-08-2011-17h17/gapless-pt1.mat');
% Vrf =  data.EL_RF;
% 
% [Irf Jrf Krf] = findsaddle(Vrf,data.X,data.Y,data.Z,2,0);
% X = (data.X); Y = (data.Y); Z = (data.Z);
% [y x z] = meshgrid(Y,X,Z);
% [Xrf Yrf Zrf] = exactsaddle(Vrf,X,Y,Z,2,0);
% 
% Qrf = spherharmxp(Vrf,Xrf,Yrf,Zrf,L,X,Y,Z);  
% VrfregenBem=spherharmcmp(Qrf,Xrf,Yrf,Zrf,L,Xregen,Yregen,Zregen);
% 
% [Ex,Ey,Ez] = gradient(VrfregenBem,1e-3*regenspacing,1e-3*regenspacing,1e-3*regenspacing);
% %[Ex,Ey,Ez] = gradient(Vrf,1e-3*data.grid(4),1e-3*data.grid(5),1e-3*data.grid(6));
% Esq0bem = Ex.^2 + Ey.^2 + Ez.^2;
% esqattrapcenter_height0 = e*Vamp^2*Esq0bem(pos,:,pos)/(4*40*mp*Omega^2)
% 
% Vrfbem = Vrf;
% 
% 
% 
% 
% 
% 
% 
% load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-08-2011-16h58/Mike_gapless_grid-pt1.mat');
% Vrf =  data.EL_RF;
% 
% [Irf Jrf Krf] = findsaddle(Vrf,data.X,data.Y,data.Z,2,0);
% X = (data.X); Y = (data.Y); Z = (data.Z);
% [y x z] = meshgrid(Y,X,Z);
% [Xrf Yrf Zrf] = exactsaddle(Vrf,X,Y,Z,2,0);
% 
% Qrf = spherharmxp(Vrf,Xrf,Yrf,Zrf,L,X,Y,Z);  
% 
% VrfregenMike=spherharmcmp(Qrf,Xrf,Yrf,Zrf,L,Xregen,Yregen,Zregen);
% 
% [Ex,Ey,Ez] = gradient(VrfregenMike,1e-3*regenspacing,1e-3*regenspacing,1e-3*regenspacing);
% %[Ex,Ey,Ez] = gradient(Vrf,1e-3*data.grid(4),1e-3*data.grid(5),1e-3*data.grid(6));
% Esq0mike = Ex.^2 + Ey.^2 + Ez.^2;
% esqattrapcenter_height0_mike = e*Vamp^2*Esq0mike(pos,:,pos)/(4*40*mp*Omega^2)
% 
% Vrfmike = Vrf;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% h = figure;
% axes1 = axes('Parent',h,'FontSize',12);
% box(axes1,'on');
% hold on
% xlabel('y (mm)');
% ylabel('U_{ps} (eV)');
% %xlim([0.14 0.45])
% %ylim([0 1*1e-6])
% plot(Yregen, esqattrapcenter_height0,'k');
% plot(Yregen, esqattrapcenter_height0_mike,'k');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% U = e^2*Vamp^2*Esq0bem/(4*40*mp*Omega^2);
% 
% [Irf Jrf Krf] = findsaddle(VrfregenBem,Xregen,Yregen,Zregen,2,0);
% X = Xregen; Y = Yregen; Z = Zregen;
% [y x z] = meshgrid(Y,X,Z);
% [Xrf Yrf Zrf] = exactsaddle(VrfregenBem,Xregen,Yregen,Zregen,2,0)
% 
% Uxy = U(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf); 
% MU = max(max(Uxy)); dL = (x(Irf+3,Jrf,Krf)-x(Irf,Jrf,Krf));
% Uxy = Uxy/MU;
% xr = (x(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-x(Irf,Jrf,Krf))/dL ;
% yr = (y(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-y(Irf,Jrf,Krf))/dL;
% [C1 C2 theta] = p2d(Uxy,xr,yr);                        % look at TestMyFunctions for testiing
% 
% fybem = (1e3/dL)*sqrt(2*C2*MU/(40*mp))/(2*pi)
% 
% 
% 
% 
% 
% 
% 
% U = e^2*Vamp^2*Esq0mike/(4*40*mp*Omega^2);
% 
% [Irf Jrf Krf] = findsaddle(VrfregenMike,Xregen,Yregen,Zregen,2,0);
% X = Xregen; Y = Yregen; Z = Zregen;
% [y x z] = meshgrid(Y,X,Z);
% [Xrf Yrf Zrf] = exactsaddle(VrfregenMike,Xregen,Yregen,Zregen,2,0)
% 
% Uxy = U(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf); 
% MU = max(max(Uxy)); dL = (x(Irf+3,Jrf,Krf)-x(Irf,Jrf,Krf));
% Uxy = Uxy/MU;
% xr = (x(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-x(Irf,Jrf,Krf))/dL ;
% yr = (y(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-y(Irf,Jrf,Krf))/dL;
% [C1 C2 theta] = p2d(Uxy,xr,yr);                        % look at TestMyFunctions for testiing
% 
% fymike = (1e3/dL)*sqrt(2*C2*MU/(40*mp))/(2*pi)
% 
% 
% 
% 
% abs(fymike-fybem)/fybem*100
% abs(1.2620720387748627*1e6-fybem)/fybem*100
% abs(fymike-1.2620720387748627*1e6)/fybem*100


% hold off

%pause;

% for i=1:length(data.Y)
%     plotVec(i)=100*abs(MikeData(i,2)-esqattrapcenter_height0(i))/max(MikeData(i,2),esqattrapcenter_height0(i));
% end
% 
% plot(data.Y,plotVec)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THE FREQUENCY OF BEMSOLVER WITHOUT SMOOTHENING THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-09-2011-10h39/gapless-pt1.mat');

% 1 mu, 180-200
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-08-2011-17h17/gapless-pt1.mat');


% 1 mu, 5 mu gap!!
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-00h54/gapless-pt1.mat');

% 1 mu, 10 mu gap!!
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-10-2011-14h37/gapless-pt1.mat');


% 1 mu, 20 mu gap!!
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-11-2011-13h07/gapless-pt1.mat');

% 1 mu, 30 mu gap!!
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-11-2011-14h10/gapless-pt1.mat');



% 1 mu, world.refine(25), 7000u ground plane, 20 mu gap!!
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-22-2011-20h08/gapless-pt1.mat');

% 1 mu, world.refine(25), 7000u ground plane, 30 mu gap!!
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-22-2011-20h29/gapless-pt1.mat');




%rotation of ground plane 15 mu gap:
% 1 mu, 15 mu gap, -5 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-14h22/gapless-pt1.mat');
% 1 mu, 15 mu gap, -4 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-15h24/gapless-pt1.mat');
% 1 mu, 15 mu gap, -2 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-13h57/gapless-pt1.mat');
% 1 mu, 15 mu gap, 0 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-11h19/gapless-pt1.mat');
% 1 mu, 15 mu gap, 0.5 deg rotation
    load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-16h21/gapless-pt1.mat');
% 1 mu, 15 mu gap, 1 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-16h03/gapless-pt1.mat');
% 1 mu, 15 mu gap, 2 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-13h36/gapless-pt1.mat');
% 1 mu, 15 mu gap, 3 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-15h43/gapless-pt1.mat');
% 1 mu, 15 mu gap, 5 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-13h02/gapless-pt1.mat');
% 1 mu, 15 mu gap, 4 deg rotation
    %load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-15h05/gapless-pt1.mat');

    
    
%rotation of ground plane 10 mu gap:    
% 1 mu, 10 mu gap!!
    load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-10-2011-14h37/gapless-pt1.mat');
% 1 mu, 10 mu gap!! 1 deg
    load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-26-2011-11h52/gapless-pt1.mat');
% 1 mu, 10 mu gap!! 2 deg
    load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-26-2011-12h14/gapless-pt1.mat');
% 1 mu, 10 mu gap!! 3 deg
    load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-26-2011-12h44/gapless-pt1.mat');
% 1 mu, 10 mu gap!! 4 deg
    load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-26-2011-14h14/gapless-pt1.mat');
% 1 mu, 10 mu gap!! 5 deg
    load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-26-2011-14h35/gapless-pt1.mat');

    
    
% 1 mu, 5 mu gap!!
load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-00h54/gapless-pt1.mat');
% 1 mu, 5 mu gap 1 deg
load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-28-2011-11h40/gapless-pt1.mat');
% 1 mu, 5 mu gap 2 deg
load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-28-2011-11h58/gapless-pt1.mat');
% 1 mu, 5 mu gap 3 deg
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-28-2011-h/gapless-pt1.mat');



Vrf =  data.EL_RF;
VrfBem = Vrf;

[Ex,Ey,Ez] = gradient(Vrf,1e-3*data.grid(4),1e-3*data.grid(5),1e-3*data.grid(6));
Esq0bemNotsmooth = Ex.^2 + Ey.^2 + Ez.^2;

U = e^2*Vamp^2*Esq0bemNotsmooth/(4*40*mp*Omega^2);


[Irf Jrf Krf] = findsaddle(Vrf,data.X,data.Y,data.Z,2,0);
X = data.X; Y = data.Y; Z = data.Z;
[y x z] = meshgrid(Y,X,Z);
[Xrf Yrf Zrf] = exactsaddle(Vrf,X,Y,Z,2,0);

Uxy = U(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf); 
MU = max(max(Uxy)); dL = (x(Irf+3,Jrf,Krf)-x(Irf,Jrf,Krf));
Uxy = Uxy/MU;
xr = (x(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-x(Irf,Jrf,Krf))/dL ;
yr = (y(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-y(Irf,Jrf,Krf))/dL;
[C1 C2 theta] = p2d(Uxy,xr,yr);                        % look at TestMyFunctions for testiing

fybem = (1e3/dL)*sqrt(2*C2*MU/(40*mp))/(2*pi);

%abs(fybem_notsmooth-1.2620720387748627*1e6)/fybem_notsmooth*100;

pos=round(size(data.Y,1)/2);
[I C] = min(U(pos,:,pos));
data.Y(C)
UBem=U(pos,:,pos);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE THE FREQUENCY OF MATHEMATICA WITHOUT SMOOTHENING THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0.2 mu,180-200
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-08-2011-18h49/Mike_gapless_grid-pt1.mat');

% 0.25 mu,185-195
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-09-2011-11h07/Mike_gapless_grid-pt1.mat');

% 1 mu, 180-200
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-09-2011-16h15/Mike_gapless_grid-pt1.mat');

% 1 mu, 180-200, 10 mu Gap, interpolation n=1
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-10-2011-18h30/Mike_gapless_grid-1-pt1.mat');

% 1 mu, 180-200, 10 mu Gap, interpolation n=5
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-10-2011-18h33/Mike_gapless_grid-5-pt1.mat');

% 1 mu, 180-200, 10 mu Gap, interpolation n=10
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-10-2011-18h34/Mike_gapless_grid-10-pt1.mat');

% 1 mu, 180-200, 10 mu Gap, interpolation n=2
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-10-2011-18h39/Mike_gapless_grid-2-pt1.mat');

% 1 mu, 180-200, 10 mu Gap, interpolation n=50
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-10-2011-20h39/Mike_gapless_grid-50-pt1.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1 mu, 5 mu Gap, interpolation n=20
load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-28-2011-11h46/Mike_gapless_grid-20-gap5-pt1.mat');

% 1 mu, 10 mu Gap, interpolation n=20
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-10-2011-18h35/Mike_gapless_grid-20-pt1.mat');

% 1 mu, 15 mu Gap, interpolation n=20
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-23-2011-11h05/Mike_gapless_grid-20-gap15-pt1.mat');

% 1 mu, 20 mu Gap, interpolation n=20
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-11-2011-12h00/Mike_gapless_grid-20-gap20-pt1.mat');

% 1 mu, 30 mu Gap, interpolation n=20
%load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/03-11-2011-14h10/Mike_gapless_grid-20-gap30-pt1.mat');



Vrf =  data.EL_RF;
VrfMike = Vrf;

[Ex,Ey,Ez] = gradient(Vrf,1e-3*data.grid(4),1e-3*data.grid(5),1e-3*data.grid(6));
Esq0MikeNotsmooth = Ex.^2 + Ey.^2 + Ez.^2;

U = e^2*Vamp^2*Esq0MikeNotsmooth/(4*40*mp*Omega^2);

[Irf Jrf Krf] = findsaddle(Vrf,data.X,data.Y,data.Z,2,0);
X = data.X; Y = data.Y; Z = data.Z;
[y x z] = meshgrid(Y,X,Z);
[Xrf Yrf Zrf] = exactsaddle(Vrf,X,Y,Z,2,0);

Uxy = U(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf); 
MU = max(max(Uxy)); dL = (x(Irf+3,Jrf,Krf)-x(Irf,Jrf,Krf));
Uxy = Uxy/MU;
xr = (x(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-x(Irf,Jrf,Krf))/dL ;
yr = (y(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-y(Irf,Jrf,Krf))/dL;
[C1 C2 theta] = p2d(Uxy,xr,yr);                        % look at TestMyFunctions for testiing

fyMike = (1e3/dL)*sqrt(2*C2*MU/(40*mp))/(2*pi);

%abs(fyMike_notsmooth-1.2620720387748627*1e6)/fyMike_notsmooth*100;
pos=round(size(data.Y,1)/2);
[I C] = min(U(pos,:,pos));
data.Y(C)
UMike=U(pos,:,pos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE THE POTENTIALS OF BEM VS. ANALYTICAL SOLTION DIRECTLY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pos=round(size(data.Y,1)/2);
% 
% aa=Vrfbem(pos,:,pos);
% bb=Vrfmike(pos,:,pos);
% 
% plot(abs(aa-bb)./aa*100)

plot(data.Y,UBem/e)
hold all
plot(data.Y,UMike/e)

[fybem fyMike];

hold off


h = figure;
plot(data.Y,VrfBem(pos,:,pos))
hold all
plot(data.Y,VrfMike(pos,:,pos))

mean(abs(VrfMike(pos,:,pos)-VrfBem(pos,:,pos))./(min(VrfMike(pos,:,pos),VrfBem(pos,:,pos))))

mean(abs(UMike/e-UBem/e));
hold off

%plot(VrfBem(pos,:,pos)./VrfMike(pos,:,pos))



% analyticPlotData=importdata('Mike_gapless_plot.dat');
% 
% h = figure;
% axes1 = axes('Parent',h,'FontSize',12);
% box(axes1,'on');
% hold(axes1,'all');
% hold all;
% xlabel({'y (\mum)'});
% ylabel({'U_{ps} (eV)'});
% xlim([140 420])
% ylim([min(analyticPlotData(:,2)) max(analyticPlotData(:,2))])
% plot(analyticPlotData(:,1),analyticPlotData(:,2),'k')
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0.05 5.12 3.84])
% set(gcf, 'PaperSize', [4.8 3.7]);
% print(h,'-dpdf','/Users/Gebhard/Documents/Berkeley/MATLAB/post-processed-simulation-data/gapless/AnalyticSol_rftrap.pdf')
% hold off; 


%error bem vs interpolated analytical approach for different number of segments
h1=figure;   
axes1 = axes('FontSize',12);
    box(axes1,'on');
    %hold(axes1,'all');
   hold all;
errorPlotX=[1 2 5 10 20 50];
errorPlotData=[6.4561e-05 1.4570e-04 1.9438e-04 2.1061e-04 2.1872e-04 2.2359e-04];
plot(errorPlotX,errorPlotData)
xlabel('Number of segments, d=10 \mum')
ylabel('Error')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0.05 5.12 3.84])
set(gcf, 'PaperSize', [4.8 4]);
hold off
print('-dpdf','error_10_mu_bemvsanalytic_segments.pdf')



% error bem vs interpolatated gapless for different gap sizes at n = 20
% world.refine(30)
h1=figure;   
axes1 = axes('FontSize',12);
    box(axes1,'on');
    %hold(axes1,'all');
   hold all;
errorPlotX=[0 5 10 15 20 30];
errorPlotData=[9.2459e-05 0.0015 2.1872e-04 5.2148e-04 1.3e-3  3.0934e-04];
plot(errorPlotX,errorPlotData)
xlabel('Gap size, n=20')
ylabel('Error')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0.05 5.12 3.84])
set(gcf, 'PaperSize', [4.8 4]);
hold off
print('-dpdf','error_bemvsanalytic_gaps.pdf')



%error bem vs interpolatated gapless for different gap sizes at n = 20
% %world.refine(25)
% h1=figure;   
% axes1 = axes('FontSize',12);
%     box(axes1,'on');
%     %hold(axes1,'all');
%    hold all;
% errorPlotX=[0 5 10 20 30];
% errorPlotData=[9.2459e-05 7.3493e-05 1.8039e-05  9.5157e-06 1.3954e-04];
% plot(errorPlotX,errorPlotData)
% xlabel('Gap size, n=20')
% ylabel('Error')
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0.05 5.12 3.84])
% set(gcf, 'PaperSize', [4.8 4]);
% hold off
% print('-dpdf','error_bemvsanalytic_gaps.pdf')

gaps=[0 5 10 15 20 30];
freqMike=[1.2563 1.2954 1.3335 1.3702 1.4052 1.4685];
freqBem=[1.2633 1.2848 1.3400 1.3613 1.3800 1.4833];

plot(gaps,freqMike)
hold all;
plot(gaps,freqBem)


%plot(abs(freqMike-freqBem)./min(freqMike,freqBem))


%error bem vs. anal. for different rotations for 15 mu gap.
rot=[-5 -4 -2 0 0.5 1 2 3 4 5];
err=[0.001 8.8079e-04 5.0878e-04 5.2148e-4 7.4587e-04 8.7643e-04 5.0798e-04 4.2077e-04 8.8079e-04 0.001];
meanerr15mu=mean(err);
meanerr=ones(size(rot,2),1)*mean(err);
figure;
plot(rot,err)
hold all;
plot(rot,meanerr)


%10 mu
rot=[0 1 2 3 4 5];
err=[1.9372e-04  7.2871e-04 7.8654e-04 8.1511e-04 7.2904e-04 7.3526e-04];
meanerr10mu=mean(err);
meanerr=ones(size(rot,2),1)*mean(err);
figure;
plot(rot,err)
hold all;
plot(rot,meanerr)

%5 mu
rot=[0 1 2];
err=[0.0015 0.0016 0.0013];
meanerr5mu=mean(err);
meanerr=ones(size(rot,2),1)*mean(err);
figure;
plot(rot,err)
hold all;
plot(rot,meanerr)



rot=[5 10 15];
err=[meanerr5mu meanerr10mu meanerr15mu];
figure;
plot(rot,err)
hold all;




