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

% LOCAL BEMSOLVER
datapathin = '/Users/Gebhard/Documents/Berkeley/Bemsolver/Sandia/bemsolver-out/';
% NETWORK BEMSOLVER LAB-25 PC
datapathin = '/Volumes/data_share/cct/Gebhard/Bemsolver-Win7/';

datapathout = '/Users/Gebhard/Documents/Berkeley/MATLAB/data/sandia/';

timenow = datestr(now,'mm-dd-yyyy-HHhMM');

mkdir(sprintf('%s%s',datapathout,timenow));

datapathout = [sprintf('%s%s',datapathout,timenow) '/'];

% specify name of txt files (the number will be added at the end!)
%dataNames = '20_3d_grid_freq_pt-400-400-21points-65refine-pt';
dataNames = '20_3d_grid_freq_pt-640-640-21points-65refine-pt';

%dataNames = '20_3d_grid_freq_pt-640-640-21points-65refine-pt'

dataNames = '20_3d_grid_freq_2mu_pt-640-640-21points-65refine-pt'



dataNames = 'gapless-pt'
%datapathin = '/Users/Gebhard/Documents/Berkeley/MATLAB/data/';
%dataNames = 'Mike_gapless_grid-20-gap5-pt'





NUM_ELECTRODES = 1; % number of non-ground electrodes
NUM_AXIS = 21; % specify number of data points per axis

NUM_DC = 40; 
NUM_Center = 2;

%UNITS: if drawing (kilian) dimensions are in microns and output in mm then
%m = 1000;
%if drawing is in mm, m =1;
m = 1000;

%COORDINATES Nikos code uses y- height, z - axial, x - radial
%if drawing uses x - axial, y - radial, z - height, use perm = [2 3 1]
%if drawing uses Nikos convention, use perm = [1 2 3]
perm = [2 3 1];
%perm = [1 2 3];


% What is the total of overlapping data structures?

nStart = 1;
nMatTot = 1;

for iterationNumber=nStart:nMatTot

    fprintf(['Importing ',sprintf('%s%i',dataNames,iterationNumber),'.txt \n']);
    
    data = load([sprintf('%s',datapathin),sprintf('%s',dataNames), sprintf('%d',iterationNumber), '.txt']); 

    %find the x y and z grids
    y = 0;
    x = 0;
    z = data(1:NUM_AXIS, 3)/m;
    for i=0:(NUM_AXIS-1)
        y(i+1) = data(NUM_AXIS*i+1, 2)/m;
        x(i+1) = data(NUM_AXIS^2*i +1, 1)/m;
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
    deltax = (xmax - xmin)/(NUM_AXIS-1);
    deltay = (ymax - ymin)/(NUM_AXIS-1);
    deltaz = (zmax - zmin)/(NUM_AXIS-1);

    %loads all the voltages and E vector into struct using dynamic naming
    for el=0:(NUM_ELECTRODES-1)
         for i=0:(NUM_AXIS-1)
            for j = 0:(NUM_AXIS-1)
                lb = NUM_AXIS^3*el + NUM_AXIS^2*i + NUM_AXIS *j + 1; %lower bound
                ub = NUM_AXIS^3*el + NUM_AXIS^2*i + NUM_AXIS *j +NUM_AXIS; %upper bound .
                struct.(['EL_phi' num2str(el)])(i+1,j+1,1:NUM_AXIS)=data(lb:ub, 4);
                
                % if loop by Gebhard, Oct 2010
                if (size(data,2)>4)
                    % i.e. Ex,Ey,Ez are calculated in bemsolver (old
                    % version), fast
                    struct.(['EL_Ex' num2str(el)])(i+1,j+1,1:NUM_AXIS)=data(lb:ub, 5);
                    struct.(['EL_Ey' num2str(el)])(i+1,j+1,1:NUM_AXIS)=data(lb:ub, 6);
                    struct.(['EL_Ez' num2str(el)])(i+1,j+1,1:NUM_AXIS)=data(lb:ub, 7);
                else
                    % i.e. Ex, Ey, Ez are NOT calculated in bemsolver (slow
                    % bemsolver, more exact). Erf will be calculated by the
                    % numerical gradient in ppt2.m
                    struct.(['EL_Ex' num2str(el)])(i+1,j+1,1:NUM_AXIS)=0;
                    struct.(['EL_Ey' num2str(el)])(i+1,j+1,1:NUM_AXIS)=0;
                    struct.(['EL_Ez' num2str(el)])(i+1,j+1,1:NUM_AXIS)=0;
                end
                
            end
        end
        struct.(['EL_phi' num2str(el)]) = permute(struct.(['EL_phi' num2str(el)]), perm);
        struct.(['EL_Ex' num2str(el)]) = permute(struct.(['EL_Ex' num2str(el)]), perm);
        struct.(['EL_Ey' num2str(el)]) = permute(struct.(['EL_Ey' num2str(el)]), perm);
        struct.(['EL_Ez' num2str(el)]) = permute(struct.(['EL_Ez' num2str(el)]), perm);

    end

    clear data
    %hard-codes to reorganize struct for specific electtode configuration
    data.project = 'design3';
    data.date_started = datestr(now);

    data.X = x;
    data.Y = y;
    data.Z = z;
    data.grid = [xmin ymin zmin deltax deltay deltaz];
    data.date_done = datestr(now);
    
 
    data.EL_RF = struct.EL_phi0;
    
   

    save([sprintf('%s',datapathout), sprintf('%s',dataNames), sprintf('%d',iterationNumber), '.mat'],'data');


end

