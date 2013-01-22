% anlz_mlp
% Generic trap analysis template. Coordianates more specific tools.
% Open up the data of a simulated trap file, get the voltages qfp[ uses for 
% a given set of multipoles and call ppt2 to do the rest.
% data is a structure with trap simulation information (potentials etc.) 
% params is a structure with all the trap parameters: voltages, RF setup, 
% ion position, stray pseudofield, frequencies, axes tilt, depth, quadrupole
% coefficients, alpha and q parameters
% main sequence is:
% run getthedata
% run regenthedata
% run setdc
% run ppt2
% run printtofile2
%
% Nikos July 2009-2010
% manipulated by Gebhard Oct 2010
%
% NOTICE: to convert bemsolver txt's into matlab data structure, run
% importd.m first!


clear all

%SIunits;

%% Inputs
%Files & Paths

% % A trap without fingers 
% importDate = 'A_trap_el4toel7-pt_03-20-2012-15h20'; % Imported .MAT file path
% dataNames = 'A_trap_el4toel7-pt';                   % Imported .MAT file name
% zMin = 1055; % FREQ   2110 for A trap
% zMax = 1795; % FREQ   3570 for A trap
% zStep = 20; % Length covered by each .MAT files

%A trap with fingers 
importDate = 'D_wire-50_100-120height_field_01-11-2013-16h48'; % Imported .MAT file path
dataNames = importDate(1:length(importDate)-17); 
zMin = -5; % FREQ
zMax = 5; % FREQ
zStep = 10; % Length covered by each .MAT files




datapath = sprintf('/home/soenke/Documents/Mathlab/trap_simulations/results/%s/',importDate);
outpath = '/home/soenke/Documents/Mathlab/trap_simulations/post-processed-simulation-data/';  %Location to save the output of this program
newfilename = '05-17-2011-13h08';   %Folder to save the output of thie program
% Trap discription 
NUM_DC = 24;            % Number of DC electrodes
NUM_Center = 0;         % Number of Center electrodes
NUM_Finger_regions = 0; % This should be 0 or equal to NUM_Center

% To chose between all mltipoles(1) and only U2(0) 
allmultipoles = 1;

% To add any kind of constraints on the way the electrodes are connected 0
% for no constrain, for others look in regenthedata
truncVoltages = 0;


pos = 0;      % 5th           %Position of the Ion along the Z axis (microns)  2538 for D trap 5th, 1270 for A trap 5th
% pos = 1065;      % 4th           %Position of the Ion along the Z axis (microns)  2538 for D trap 5th, 1270 for A trap 5th
% pos = 760;      % 3rd           %Position of the Ion along the Z axis (microns)  2538 for D trap 5th, 1270 for A trap 5th
% pos = 2085;      % 8th           %Position of the Ion along the Z axis (microns)  2538 for D trap 5th, 1270 for A trap 5th


% In case of external voltage set
extV = 0; % 1 if yes (External voltage set 'Vext' should be added in the External voltage portion below)

%% External Voltages
if extV
    Vext = [0;-0.272400066597325;1.15638746689082;2.75468878071420;-11.5353564487889;2.76701343278494;1.07421239243774;-0.372809107750100;-0.317891444873378;0;0;0;-0.490177618040699;-0.958952250010321;-4.98473231177900;-3.75881722075925;-4.78237006933698;-0.943080796521006;-0.504297840329082;-0.321827523335782;0;0;-0.574791231420775];
end

%% Core

RF_offset = 0;

position = pos/1000;
data.pre = 0;

data = getdata(datapath,position,dataNames,1,zMin,zMax,zStep,NUM_DC,NUM_Center);
data = regenthedata(data,0,0,position,9,NUM_DC,NUM_Center,truncVoltages,RF_offset);
data.allmultipoles = allmultipoles;
data = trapknobs(data,position,0,0,1);

%params.datesim =  simdate;
%params.datedata = data'Y:/cct/Sankar/MATLAB/trap simulationdate;
params.daterun = datestr(date);
params.position = position;
ind = true; 
reg = false;
az = 4.5e-3;
ax = -0.002;
phi = 0;


Ex = 0;
Ey = 0;
Ez = 0;

params.rfamplitude = 240; 
params.frequency = 60e6;
U1 = -.22;
U2 = 6.0; 
U3 = .22;  
U4 = 0; 
U5 = 0;
params.E = [Ex Ey Ez];
params.scale = 1;
if extV
    el = Vext;
else
    el = setdc(data,params.frequency,-Ex,-Ey,-Ez,U1,U2,U3,U4,U5,ax,az,phi,ind,reg);
end
params.VELDC = el;
%params.VELDC = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];

outpath2 = sprintf('%s/',outpath);
newfilename2 = 0;
params = ppt2(params,data,'N',NUM_DC,NUM_Center,2,0,0,2,outpath2,newfilename2,pos,truncVoltages,RF_offset);

mes = sprintf('The secular frequencies are: (%G, %G, %G) Hz.\n',params.f(1),params.f(2),params.f(3));
disp(mes);
mes = sprintf('The trapdepth is: %G eV \n',params.trapdepth);
disp(mes);
mes = sprintf('The ion sits at (%G,%G,%G) micron.\n',1e3*params.ionpos(1),1e3*params.ionpos(2),1e3*params.ionpos(3));
disp(mes);

save('/home/soenke/Documents/Mathlab/trap_simulations/traim.mat');

% CC = reshape(data.C,1,184);
% CC = CC';
% for a = 1:5
%     CC(:,a) = CC(:,1);
%     CC(185,a) = a;
% end
CC=data.C
dlmwrite('/home/soenke/Documents/Mathlab/trap_simulations/D_trap/Cfiles/D_trap_mid.txt', CC, 'delimiter', ' ')



