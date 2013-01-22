% just copy and paste to the workspace
% Nikos 2008-2009

% tensor*vector
T3 = zeros(3,3,3);
T3(:,:,1)= eye(3);
T = zeros(3,3,3,3);
T(:,:,:,1) = T3;
V = [0;0;1];
PR = multiply(T,V);

% findsaddle
[X Y Z] = meshgrid(-0.8:0.1:0.8,-0.9:0.1:0.9,-1:0.1:1);
grid = [-0.8 -0.9 -1.0 0.1 0.1 0.1];
f = X.^2+3*Y.^2+2*Z.^2+X.^3-Y.^3-Z.^3-5*X.^4;
tic;
[I J K] = findsaddle(f,grid,3);
toc
PlotPot(f,I,J,K,grid,2);


%exactsaddle
for i=1:100
Xax = -1:0.1:1; Yax = -1:0.1:1; Zax = -1:0.1:1;
grid = [Xax(1) Yax(1) Zax(1) Xax(2)-Xax(1) Yax(2)-Yax(1) Zax(2)-Zax(1)];
[Y X Z] = meshgrid(Yax,Xax,Zax);
f = 10*((X-0.12).^2+(Y-0.3975).^2-2*(Z-0.3).^2)+((X-0.12).^2-(Y-0.397).^2)+i*X/100; % 3d 
%f = 10*((X-0.12).^2-(Y-0.001).^2)+i*X/100; % 2d
%tic;
[I J K] = polysaddle(f,Xax,Yax,Zax,3,0);
%toc
x(i)=I; y(i)=J; z(i)=K;
end
plot(y);

%p2d rotated 2d polynomial fit
Xax = -1:0.1:1; Yax = -1:0.1:1;
[X Y] = meshgrid(Xax,Yax);
theta = 30*pi/180;
Xr = X*cos(theta)+Y*sin(theta);
Yr = -X*sin(theta)+Y*cos(theta);
F = 1*Xr.^2+0.5*Yr.^2;
[Af Bf th]=p2d(F,X,Y);

%rotated mulitpoles
for i=1:180
Xax = -1:0.1:1; Yax = -1:0.1:1; Zax = -1:0.1:1;
[Y X Z] = meshgrid(Yax,Xax,Zax);
theta = (i-90)*pi/180;
Xr = X*cos(theta)+Y*sin(theta);
Yr = -X*sin(theta)+Y*cos(theta);
Zr = Z;
F = 1*(Xr.^2-Yr.^2)/2;
Q = spherharmxp(F,0,0,0,3,Xax,Yax,Zax);
A(i) = 2*sqrt( (3*Q(8))^2+(3*Q(9))^2 );
Theta(i) = 45*(sign(Q(9)))-90*atan((3*Q(8))/(3*Q(9)))/pi;
%Theta(i) = -90*atan((3*Q(8))/(3*Q(9)))/pi;
end
plot(A); pause
plot(-89:90,(Theta));pause;

%spherharmxp
clear; 
Xax = -2:0.2:2; Yax = -2:0.2:2; Zax = -2:0.2:2;
[Y X Z] = meshgrid(Yax,Xax,Zax);
%F = 0.5*(2*(Z-0.05).^3-3*((X-0.05).^2+(Y-0.05).^2).*(Z-0.05));
%F = -1.5*(4*(X-0.05).*(Z-0.05).^2-(X-0.05).*((X-0.05).^2+(Y-0.05).^2));
%F = -1.5*(4*(Y-0.05).*(Z-0.05).^2-((X-0.05).^2+(Y-0.05).^2).*(Y-0.05));
%F = 15*((Z-0.05).*((X-0.05).^2-(Y-0.05).^2));
%F = 30*(X-0.05).*(Y-0.05).*(Z-0.05);
%F = 15*[-(X-0.05).^3+3*(X-0.05).*(Y-0.05).^2];
F = 15*[(Y-0.05).^3-3*(Y-0.05).*(X-0.05).^2];
%F = (Z-0.05).^2-(X-0.05).^2/2-(Y-0.05).^2/2+(X-0.05).*(Y-0.05)+[(X-0.05).^3-3*(X-0.05).*(Y-0.05).^2];
Q = spherharmxp(F,0.05,0.05,0.05,5,Xax,Yax,Zax);

% trapdepth
Xax = -1:0.05:1; Yax = -1:0.05:1; Zax = -1:0.05:1;
[Y X Z] = meshgrid(Yax,Xax,Zax);
F = Z.^2+X.^2+Y.^2-X.^3;
[D Xe Ye Ze] = trapdepth(F,Xax,Yax,Zax,11,11,11);

% getthedata
datapath = 'D:\User\x22473\Nikos\Matlab\cpo-trap-simulations\cpo_data';
simdate = '090408';
z = 0.39:0.01:0.62;
for i = 1:size(z,2)
    data = getthedata_ref(datapath,simdate,z(i));
    Zmin(i) = min(data.Z);
    Zmax(i) = max(data.Z);
    Length(i) = size(data.Z,1);
end
plot(Zmin); pause; plot(Zmax); pause; plot(Length); 

% mplcoeff and trapknobs
data.X = -1.5:0.1:1.5; data.Y = -1.5:0.1:1.5; data.Z = -1.5:0.1:1.5;
data.grid = [-1.5 -1.5 -1.5 0.1 0.1 0.1];
[Y X Z] = meshgrid(data.Y,data.X,data.Z);
data.Vrf = (X-0.5).^2/2-(Y+0.2).^2/2;
data.W1 = X-0.5; 
data.W2 = Y+0.5; 
data.W3 = Z-0.1; 
data.W4 = (X-0.5).^2/2-(Y+0.5).^2/2;
data.W5 = (Z-0.1).^2-(X-0.5).^2/2-(Y+0.5).^2/2; 
data.W6 = (X-0.5).*(Y+0.5)/2;
data.W7 = (Y+0.5).*(Z-0.1)/2; 
data.W8 = (Z-0.1).*(X-0.5)/2; 
data.W9 = 2*(X-0.5); 
data.W10 = 2*(Y+0.5);
data.N1 = 2*(Z-0.1); 
data.N2 = (X-0.5).^2-(Y+0.5).^2; 
data.N3 = 2*(Z-0.1).^2-(X-0.5).^2-(Y+0.5).^2; 
data.N4 = (X-0.5).*(Y+0.5);
data.N5 = (Y+0.5).*(Z-0.1); 
data.N6 = (Z-0.1).*(X-0.5); 
data.N7 = 3*(X-0.5); 
data.N8 = 3*(Y+0.5); 
data.N9 = 3*(Z-0.1); 
data.N10 = 3*(X-0.5).^2/2-3*(Y+0.5).^2/2; 
data.Vc = 3*(Z-0.1).^2-3*(X-0.5).^2/2-3*(Y+0.5).^2/2;
position = 0.1;
M = mplcoeff(data,position,8);                                    % find multipole coefficients for all electrodes
rs = 1; s = [1];                                                % rescale the coefficients from mm to some lengthscale eg 1 micron <-> rs = 1e-3
for i=1:8
    for j = 1:2*i+1
        s = horzcat(s,rs^i);
    end
end
S = diag(s);
M = S*M;
mesh(M);
data.M = M;
data = trapknobs(data,position);

% quadparams
data.X = -1.5:0.1:1.5; data.Y = -1.5:0.1:1.5; data.Z = -1.5:0.1:1.5;
[Y X Z] = meshgrid(data.Y,data.X,data.Z);
%U = (X-0.5).^2/2-(Y+0.5).^2/2;
%U =  (Z-0.1).^2-(X-0.5).^2/2-(Y+0.5).^2/2;
%U = (X-0.5).*(Y+0.5)/2;
U =  (Y+0.5).*(Z-0.1)/2; 
%U =  (Z-0.1).*(X-0.5)/2; 
quadparams(U,data,1);

% ppt2
data.grid = [-1.5 -1.5 -1.5 0.1 0.1 0.1];
data.X = (-1.5:0.1:1.5)'; data.Y = (-1.5:0.1:1.5)'; data.Z = (-1.5:0.1:1.5)';
[Y X Z] = meshgrid(data.Y,data.X,data.Z);
data.Vrf = (X-0.5).^2/2-(Y+0.5).^2/2+(X-0.5).*(Y+0.5)/2;
data.Erf = sqrt((X-0.5).^2+(Y+0.5).^2);
data.W1 = 0;%X-0.5; 
data.W2 = 0;%Y+0.5; 
data.W3 = 0;%Z-0.1; 
data.W4 = (X-0.5).^2/2-(Y+0.5).^2/2;
data.W5 = (Z-0.1).^2-(X-0.5).^2/2-(Y+0.5).^2/2; 
data.W6 = (X-0.5).*(Y+0.5)/2;
data.W7 = 0;%(Y+0.5).*(Z-0.1)/2; 
data.W8 = 0;%(Z-0.1).*(X-0.5)/2; 
data.W9 = 0;%2*(X-0.5); 
data.W10 = 0;%2*(Y+0.5);
data.N1 = 0;%2*(Z-0.1); 
data.N2 = 0;%(X-0.5).^2-(Y+0.5).^2; 
data.N3 = 0;%2*(Z-0.1).^2-(X-0.5).^2-(Y+0.5).^2; 
data.N4 = 0;%(X-0.5).*(Y+0.5);
data.N5 = 0;%(Y+0.5).*(Z-0.1); 
data.N6 = 0;%(Z-0.1).*(X-0.5); 
data.N7 = 0;%3*(X-0.5); 
data.N8 = 0;%3*(Y+0.5); 
data.N9 = 0;%3*(Z-0.1); 
data.N10 = 0;%3*(X-0.5).^2/2-3*(Y+0.5).^2/2; 
data.Vc = 0;%3*(Z-0.1).^2-3*(X-0.5).^2/2-3*(Y+0.5).^2/2;
position = 0.1;
params.E = [0 0 0];
params.position = 0.1;
params.Ver = 0;
params.Hor = 0; 
params.W = 0.1*[0 0 0 0 1 0.2 0 0 0 0];
params.N = 0.01*[0 0 0 0 0 0 0 0 0 0];
params.Center = 0;
params.scale = 1;
params.rfamplitude = 100;
params.frequency = 15e6;
params.datesim = 'today';
params.datedata = 'today';
params = ppt2(params,data,'N');
printtofile2(1,'ind',params);