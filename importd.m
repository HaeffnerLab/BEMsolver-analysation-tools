% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% importd.m
%
% original by Mike, modified by Gebhard Oct 2010
%
% Imports bemsolver data according to Gebhard's file saving notation. Takes
% .txt files as input and converts it into a data structure .dat. Saves
% mat-structure under the same filename as the .txt.
%
% importd.m can import multiple .txt files (see for loop). The number of
% files to be imported can be adjusted via nMatTot.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
set = 1;

if (set == 1)
    datapathin = '/home/soenke/Documents/Mathlab/trap_simulations/results/'; % specify name of txt files (the number will be added at the end!)
    dataNames = 'A_wire-50x250z_89-139_WR65_field';
elseif(set == 2)
     datapathin = '/home/soenke/Mathlab/trap_simulations/';
     dataNames = 'Eu_Trap-pt1_pos_E7E8_cor2';
elseif(set == 3)
    datapathin = '/home/soenke/Mathlab/trap_simulations/'; % specify name of txt files (the number will be added at the end!)
    dataNames = 'Eu_Trap-pt1_pos_E7E8_cor_2_2';
 end

datapathout = '/home/soenke/Documents/Mathlab/trap_simulations/results/';
timenow = datestr(now,'mm-dd-yyyy-HHhMM');
fl_name = sprintf('%s%c%s',dataNames,'_',timenow);
mkdir(sprintf('%s%s',datapathout,fl_name));
datapathout = [sprintf('%s%s',datapathout,fl_name) '/'];

data.NUM_AXIS = 31; % specify number of data points per axis
NUM_AXIS=data.NUM_AXIS;
data.NUM_ELECTRODES = 24; % number of non-ground electrodes
data.NUM_DC = 23; 
data.NUM_CENTER = 0;

%UNITS: if drawing (kilian) dimensions are in microns and output in mm then
%m = 1000 if drawing is in mm, m =1;
m = 1000;

%COORDINATES Nikos code uses y- height, z - axial, x - radial
%if drawing uses x - axial, y - radial, z - height, use perm = [2 3 1] (Euro trap)
%if drawing uses y - axial, x - radial, z - height, use perm = [1 3 2] (Sqip D trap, GG trap)
%if drawing uses Nikos convention, use perm = [1 2 3]

perm = [1 3 2];

% What is the total of overlapping data structures?

nStart = 1;
nMatTot = 1;

for iterationNumber=nStart:nMatTot

    fprintf(['Importing ',sprintf('%s%i',dataNames,iterationNumber),'.txt \n']);
    
    imdata = load([sprintf('%s',datapathin),sprintf('%s',dataNames), sprintf('%d',iterationNumber), '.txt']); 

    %find the x y and z grids
    y = 0;
    x = 0;
    z = imdata(1:data.NUM_AXIS, 3)/m;
    for i=0:(data.NUM_AXIS-1)
        y(i+1) = imdata(data.NUM_AXIS*i+1, 2)/m;
        x(i+1) = imdata(data.NUM_AXIS^2*i +1, 1)/m;
    end
    x = x';
    y = y';
    coord = [x y z];
    x = coord(:,perm(1));
    y = coord(:,perm(2));
    z = coord(:,perm(3));

    %find  min, max and spacing of the axes
    xmin = min(x);
    xmax = max(x);
    ymin = min(y);
    ymax = max(y);
    zmin = min(z);   
    zmax = max(z);
    deltax = (xmax - xmin)/(data.NUM_AXIS-1);
    deltay = (ymax - ymin)/(data.NUM_AXIS-1);
    deltaz = (zmax - zmin)/(data.NUM_AXIS-1);

    %loads all the voltages and E vector into struct using dynamic naming
    for el=0:(data.NUM_ELECTRODES-1)
         for i=0:(data.NUM_AXIS-1)
            for j = 0:(data.NUM_AXIS-1)
                lb = data.NUM_AXIS^3*el + data.NUM_AXIS^2*i + data.NUM_AXIS *j + 1; %lower bound
                ub = data.NUM_AXIS^3*el + data.NUM_AXIS^2*i + data.NUM_AXIS *j +data.NUM_AXIS; %upper bound .
                struct.(['EL_phi' num2str(el)])(i+1,j+1,1:data.NUM_AXIS)=imdata(lb:ub, 4);
                
                % if loop by Gebhard, Oct 2010
                if (size(imdata,2)>4)
                    % i.e. Ex,Ey,Ez are calculated in bemsolver (old
                    % version), fast
                    struct.(['EL_Ex' num2str(el)])(i+1,j+1,1:data.NUM_AXIS)=imdata(lb:ub, 5);
                    struct.(['EL_Ey' num2str(el)])(i+1,j+1,1:data.NUM_AXIS)=imdata(lb:ub, 6);
                    struct.(['EL_Ez' num2str(el)])(i+1,j+1,1:data.NUM_AXIS)=imdata(lb:ub, 7);
                else
                    % i.e. Ex, Ey, Ez are NOT calculated in bemsolver (slow
                    % bemsolver, more exact). Erf will be calculated by the
                    % numerical gradient in ppt2.m
                    struct.(['EL_Ex' num2str(el)])(i+1,j+1,1:data.NUM_AXIS)=0;
                    struct.(['EL_Ey' num2str(el)])(i+1,j+1,1:data.NUM_AXIS)=0;
                    struct.(['EL_Ez' num2str(el)])(i+1,j+1,1:data.NUM_AXIS)=0;
                end
                
            end
        end
        struct.(['EL_phi' num2str(el)]) = permute(struct.(['EL_phi' num2str(el)]), perm);
        struct.(['EL_Ex' num2str(el)]) = permute(struct.(['EL_Ex' num2str(el)]), perm);
        struct.(['EL_Ey' num2str(el)]) = permute(struct.(['EL_Ey' num2str(el)]), perm);
        struct.(['EL_Ez' num2str(el)]) = permute(struct.(['EL_Ez' num2str(el)]), perm);

    end

    clear imdata
    %hard-codes to reorganize struct for specific electrode configuration
    data.project = 'design3';
    data.date_started = datestr(now);

    data.X = x;
    data.Y = y;
    data.Z = z;
    data.grid = [xmin ymin zmin deltax deltay deltaz];
    data.date_done = datestr(now);
    
    data.EL_RF = struct.EL_phi0;%+struct.EL_phi1;
    
    
    % Need to edit this if you are using out of phase! layer1=DC1 ect if
    % layer 2 is DC1 if +1
    for iii=1:(data.NUM_DC)
        data.(['EL_DC' num2str(iii)]) = struct.(['EL_phi' num2str(iii)]);
    end
     
    for iii=1:(data.NUM_CENTER)
        pos = data.NUM_DC+iii;
        data.(['EL_CNT' num2str(iii)]) = struct.(['EL_phi' num2str(pos)]);
    end  
    

    
%  DO NOT USE THIS WITHOUT MUCH THOUGH !!!!!!!!!!    
%     data.Erf = sqrt((struct.EL_Ex2+struct.EL_Ex3).^2 + (struct.EL_Ey2+struct.EL_Ey3).^2 + (struct.EL_Ez2+struct.EL_Ez3).^2); %sum of e fields of RF electrodes
%     data.Vrf = struct.EL_phi2 + struct.EL_phi3;
%     data.Vc1 = struct.EL_phi1; %lower DC electrode
%     data.Vc2 = struct.EL_phi2; %upper DC electrode
%     data.Vc3 = struct.EL_phi3; %
%     data.Vc4 = struct.EL_phi4; %
%     data.N1 = struct.EL_phi4;	%N us upper electrodes
%     data.N2 = struct.EL_phi5;
%     data.N3 = struct.EL_phi6;
%     data.N4 = struct.EL_phi7;
%     data.N5 = struct.EL_phi8;
%     data.N6 = struct.EL_phi9;
%     data.N7 = struct.EL_phi10;
%     data.N8 = struct.EL_phi11;
%     data.N9 = struct.EL_phi12;
%     data.N10 = struct.EL_phi13;
%     data.W1 = struct.EL_phi14;	%W is lower electrodes
%     data.W2 = struct.EL_phi15;
%     data.W3 = struct.EL_phi16;
%     data.W4 = struct.EL_phi17;
%     data.W5 = struct.EL_phi18;
%     data.W6 = struct.EL_phi19;
%     data.W7 = struct.EL_phi20;
%     data.W8 = struct.EL_phi21;
%     data.W9 = struct.EL_phi22;
%     data.W10 = struct.EL_phi23;
    

    save([sprintf('%s',datapathout), sprintf('%s',dataNames), sprintf('%d',iterationNumber), '.mat'],'data');
    %clear

    %load data.mat
    %RFanalysis
    %DCanalysis


end

Ef = data.EL_RF;
for a = 1:data.NUM_AXIS
    for b = 1:data.NUM_AXIS
        E(a,b) = Ef(a,b,data.NUM_AXIS);
    end
end
figure; surface(x,y,E)


