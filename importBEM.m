% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% importBEM.m
%
% original by Soenke April 2012
%
% Imports BEM solver *field.txt into a matlab data structure. Recognizes
% stepsize and step number
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

Num_RFEL=2; %Change this to the number of RF electrodes
Out_of_phase=1; %1=out of phase drive, 0=in phase drive



datapathin = '/home/soenke/Documents/Mathlab/trap_simulations/results/'; % specify name of txt files (the number will be added at the end!)
dataNames = 'layertrap2_field1';

datapathout = '/home/soenke/Documents/Mathlab/trap_simulations/results/';
timenow = datestr(now,'mm-dd-yyyy-HHhMM');
fl_name = sprintf('%s%c%s',dataNames,'_',timenow);
mkdir(sprintf('%s%s',datapathout,fl_name));
datapathout = [sprintf('%s%s',datapathout,fl_name) '/'];

fprintf(['Importing ',sprintf('%s%i',dataNames),'.txt \n']); 
data = load([sprintf('%s',datapathin),sprintf('%s',dataNames), '.txt']); %format datapathin+dataNames.txt

perm=[1 2 3]; % perm 123 if as in autocad, change accordingly

maxvalue=max(data(:,3));
mindex=find(data(:,3)==maxvalue);
steps=mindex(1);

deltax=data((steps^2+1),perm(1))-data(1,perm(1));
deltay=data((steps+1),perm(2))-data(1,perm(2));
deltaz=data(2,perm(3))-data(1,perm(3));

if (deltax==deltay && deltax==deltaz)
else
    fprintf('WARNING NOT a cubic simulation area!!!')
end

xmin = min(data(:,perm(1)));
xmax = max(data(:,perm(1)));
ymin = min(data(:,perm(2)));
ymax = max(data(:,perm(2)));
zmin = min(data(:,perm(3)));
zmax = max(data(:,perm(3)));

x=xmin:(xmax-xmin)/(steps-1):xmax;
y=ymin:(ymax-ymin)/(steps-1):ymax;
z=zmin:(zmax-zmin)/(steps-1):zmax;
x=x';
y=y';
z=z';

coord=[x y z];
x=coord(:,perm(1));
y=coord(:,perm(2));
z=coord(:,perm(3));

num_electrodes=length(find(data(:,3)==maxvalue))/(steps*steps);

for EL=0:(num_electrodes-1)
  
end
        


