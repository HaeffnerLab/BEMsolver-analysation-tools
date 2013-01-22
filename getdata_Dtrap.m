function data = getthedata_Dtrap(path,position,dataNames,plotOpt,zMin,zMax,zStep,NUM_DC,NUM_Center)
% generate the data structure for trap operation around axial position
% centered at "position". the consecutive data srtuctures have overlapping
% first and last points, i.e. ptI.Z(31) = ptI+1.Z(1)

% ... from microns into millimeters:
zMin = zMin/1000;
zMax = zMax/1000;
zStep = zStep/1000;


% What is the total of overlapping data structures?
% rounding even necessary if nMatTot has a integer value!!! (had a problem
% in the (I<1)||(I>nMatTot); somehow 17>17.00000 = 1 ?????
nMatTot = round((zMax-zMin)/zStep);

% define a "grid" of .mat-files
Zlim = zMin:zStep:zMax;

% find in which of the .mat-files the chosen position is centered
I = find(Zlim>position,1,'first')-1;

if (I<1)||(I>nMatTot),
    fprintf('Invalid ion position. Quitting.\n');
    return
end
N = sign(2*position-Zlim(I)-Zlim(I+1));
if (I==1)&&(N==-1),
    % Ion is in the first grid
    d = load([sprintf('%s',path),sprintf('%s',dataNames),'1.mat']); 
    data = d.data;
    return
elseif (I==nMatTot)&&(N==1),
    % Ion is in the last grid
    d = load([sprintf('%s',path),sprintf('%s',dataNames),sprintf('%i.mat',nMatTot)]); 
    data = d.data;
    return
else
    % somewhere inbetween
    d = load([sprintf('%s',path),sprintf('%s',dataNames),sprintf('%i.mat',I)]); 
    data0=d.data;
    fieldnames(data0);
    [dum K] = min(abs(data0.Z-position));
    clear data0;
    if N == 1,
        % pos lies closer to the next .mat-grid than the last grid 
        d = load([sprintf('%s',path),sprintf('%s',dataNames),sprintf('%i.mat',I)]);
        data1 = d.data;
        d = load([sprintf('%s',path),sprintf('%s',dataNames),sprintf('%i.mat',I+N)]); 
        data2 = d.data;
        [dum K] = min(abs(data1.Z-position));
        K1 = K-floor(numel(data1.Z)/2); % these are not the neatest definitions, but ok for now
        K2 = numel(data1.Z);    
        K3 = 2; 
        K4 = K-floor(numel(data1.Z)/2);
        %fprintf('K1=%i K2=%i K3=%i K4=%i\n',K1,K2,K3,K4);
    elseif N  == -1,
        % pos lies closer to the last .mat-grid than the next grid 
        d = load([sprintf('%s',path),sprintf('%s',dataNames),sprintf('%i.mat',I+N)]); 
        data1 = d.data;
        d = load([sprintf('%s',path),sprintf('%s',dataNames),sprintf('%i.mat',I)]); 
        data2 = d.data;
        K1 = K+floor(numel(data1.Z)/2);
        K2 = numel(data1.Z)-1;
        K3 = 1; 
        K4 = K+floor(numel(data1.Z)/2);
        %fprintf('K1=%i K2=%i K3=%i K4=%i\n',K1,K2,K3,K4);
    else
        d = load([sprintf('%s',path),sprintf('%s',dataNames),sprintf('%i.mat',I)]); 
        data = d.data;
        return
    end
end

data = data1;

% creates a new grid by joining grid with adjacent one (either left or right one; that is
% determined by the sign of N)

for i=K1:K2       
        
    data.Z(i-K1+1)=data1.Z(i); 
    
    data.EL_RF(:,:,i-K1+1) = data1.EL_RF(:,:,i);
    
    for iii=1:(NUM_DC)
        data.(['EL_DC' num2str(iii)])(:,:,i-K1+1) = data1.(['EL_DC' num2str(iii)])(:,:,i);
    end
     
    for iii=1:(NUM_Center)
        data.(['EL_CNT' num2str(iii)])(:,:,i-K1+1) = data1.(['EL_CNT' num2str(iii)])(:,:,i);
    end  
    
    
%     data.N1(:,:,i-K1+1) = data1.N1(:,:,i);    
%     data.N10(:,:,i-K1+1) = data1.N10(:,:,i);    
%     data.N2(:,:,i-K1+1) = data1.N2(:,:,i);
%     data.N3(:,:,i-K1+1) = data1.N3(:,:,i);
%     data.N4(:,:,i-K1+1) = data1.N4(:,:,i);
%     data.N5(:,:,i-K1+1) = data1.N5(:,:,i);   
%     data.N6(:,:,i-K1+1) = data1.N6(:,:,i);   
%     data.N7(:,:,i-K1+1) = data1.N7(:,:,i);   
%     data.N8(:,:,i-K1+1) = data1.N8(:,:,i);    
%     data.N9(:,:,i-K1+1) = data1.N9(:,:,i);
%     data.Vc(:,:,i-K1+1) = data1.Vc(:,:,i);    
%     
%     data.Erf(:,:,i-K1+1) = data1.Erf(:,:,i);
%     data.W1(:,:,i-K1+1) = data1.W1(:,:,i);   
%     data.W10(:,:,i-K1+1) = data1.W10(:,:,i);
%     data.W2(:,:,i-K1+1) = data1.W2(:,:,i);    
%     data.W3(:,:,i-K1+1) = data1.W3(:,:,i);
%     data.W4(:,:,i-K1+1) = data1.W4(:,:,i);    
%     data.W5(:,:,i-K1+1) = data1.W5(:,:,i);    
%     data.W6(:,:,i-K1+1) = data1.W6(:,:,i);
%     data.W7(:,:,i-K1+1) = data1.W7(:,:,i);    
%     data.W8(:,:,i-K1+1) = data1.W8(:,:,i);
%     data.W9(:,:,i-K1+1) = data1.W9(:,:,i);  
%     
    
    
    
end
for i=K3:K4        
    
    data.Z(i-K3+max(K2-K1,-1)+2)=data2.Z(i); 
    
    data.EL_RF(:,:,i-K3+max(K2-K1,-1)+2) = data2.EL_RF(:,:,i);
    
    for iii=1:(NUM_DC)
        data.(['EL_DC' num2str(iii)])(:,:,i-K3+max(K2-K1,-1)+2) = data2.(['EL_DC' num2str(iii)])(:,:,i);
    end
     
    for iii=1:(NUM_Center)
        data.(['EL_CNT' num2str(iii)])(:,:,i-K3+max(K2-K1,-1)+2) = data2.(['EL_CNT' num2str(iii)])(:,:,i);
    end  
    
    
%     data.Z(i-K3+max(K2-K1,-1)+2)=data2.Z(i);   
%     data.N1(:,:,i-K3+max(K2-K1,-1)+2) = data2.N1(:,:,i);    
%     data.N10(:,:,i-K3+max(K2-K1,-1)+2) = data2.N10(:,:,i);    
%     data.N2(:,:,i-K3+max(K2-K1,-1)+2) = data2.N2(:,:,i);   
%     data.N3(:,:,i-K3+max(K2-K1,-1)+2) = data2.N3(:,:,i);
%     data.N4(:,:,i-K3+max(K2-K1,-1)+2) = data2.N4(:,:,i);
%     data.N5(:,:,i-K3+max(K2-K1,-1)+2) = data2.N5(:,:,i);
%     data.N6(:,:,i-K3+max(K2-K1,-1)+2) = data2.N6(:,:,i);
%     data.N7(:,:,i-K3+max(K2-K1,-1)+2) = data2.N7(:,:,i);    
%     data.N8(:,:,i-K3+max(K2-K1,-1)+2) = data2.N8(:,:,i);
%     data.N9(:,:,i-K3+max(K2-K1,-1)+2) = data2.N9(:,:,i);
%     data.Vc(:,:,i-K3+max(K2-K1,-1)+2) = data2.Vc(:,:,i);
%     data.Vrf(:,:,i-K3+max(K2-K1,-1)+2) = data2.Vrf(:,:,i);    
%     data.Erf(:,:,i-K3+max(K2-K1,-1)+2) = data2.Erf(:,:,i);    
%     data.W1(:,:,i-K3+max(K2-K1,-1)+2) = data2.W1(:,:,i);    
%     data.W10(:,:,i-K3+max(K2-K1,-1)+2) = data2.W10(:,:,i);    
%     data.W2(:,:,i-K3+max(K2-K1,-1)+2) = data2.W2(:,:,i);
%     data.W3(:,:,i-K3+max(K2-K1,-1)+2) = data2.W3(:,:,i);
%     data.W4(:,:,i-K3+max(K2-K1,-1)+2) = data2.W4(:,:,i);    
%     data.W5(:,:,i-K3+max(K2-K1,-1)+2) = data2.W5(:,:,i);
%     data.W6(:,:,i-K3+max(K2-K1,-1)+2) = data2.W6(:,:,i);    
%     data.W7(:,:,i-K3+max(K2-K1,-1)+2) = data2.W7(:,:,i);    
%     data.W8(:,:,i-K3+max(K2-K1,-1)+2) = data2.W8(:,:,i);   
%     data.W9(:,:,i-K3+max(K2-K1,-1)+2) = data2.W9(:,:,i);  
    
    
end
data.grid = [min(data.X) min(data.Y) min(data.Z) data1.grid(4) data1.grid(5) data1.grid(6)];

    
% plot(data.Z,'--*'); 
% title('getthedata.m check: you will see a straight line if the data was generated succesfully.');
% pause;
    
    %title(sprintf('%f',position));
 
% for checking
%path = 'D:\User\x22473\Nikos\Matlab\cpo-trap-simulations\cpo_data';
%simdate = '090408';
%for i=0:1000
%    getthedata(path,simdate,4*i/1000-0.12); pause(0.3);
%end
