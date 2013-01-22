function  out = ppt2_Dtrap(params,data,func,NUM_DC,NUM_Center,dcplot,rfplot,plotPseudopotential,plotTrappotential,outpath,newfilename,posIndex,truncVoltages,RF_offset)
% out ppt2(params,data,func)
% post processing tool
% params 
%   the trap parameters. Input voltages must not include compensation
% data 
%   a cpo-simulation data-structure with the electrode potentials in the
%   trapping region
% func 
%   determines the task performed:
%       'E' determines the stray electric field for given dc voltages 
%       'C' determines compensation parameters for given starting voltages 
%           and given stray field
%       'N' do not optimize, just analyze the trap, assuming everything is
%           ok
%
% Start with input values for all the dc voltages and RF parameters. These
% must not include compensation parameters (added seperately).
% Find stray field that would be compensated in the given configuration.
% The axes order in the potential data is:
% I=radial horizontal (X); J=radial vertical (Y); K=axial (Z)
% ie I->X, J->Y, K->Z
% Nikos, January 2009


% rfplot, dcplot:   sets plotoptions...



e=1.60217646e-19;                                                                       % elementary charge unit in SI
mp=1.67262158e-27;                                                                      % proton mass in SI units
Mn = 40;                                                                                            % Mass number of the ion

qcheck = false;             % I think this means perforf quality checks
Zval = params.position;
VELDC = params.VELDC;

scale = params.scale;
out = params;

grid = data.grid;
X = normalize(data.X); Y = normalize(data.Y); Z = normalize(data.Z);
[y x z] = meshgrid(Y,X,Z);

RFampl = params.rfamplitude;                                               % RF parameters
Freq = params.frequency;                    
Omega = 2*pi*Freq;   
r0 = 1;                                                                    % lengthscale of multipole expansion in millimeters
V0 = Mn*mp*(2*pi*params.frequency)^2*(r0*1e-3)^2/e;                        % god given voltage in SI 

[Irf Jrf Krf] = findsaddle(data.EL_RF,X,Y,Z,2,Zval);
warn('RF',data.EL_RF,Irf,Jrf,Krf);
Vrf = RFampl*data.EL_RF;
plotpot(Vrf,Irf,Jrf,Krf,data.grid,rfplot,'RF potential','V_{rf} (Volt)',outpath,newfilename);

%Vdc = VDC1(scale*(W-Hor),scale*(N+Hor),scale*(Cnt+Ver),0,0,0);                                                 % DC parameters                                                            
%[Idc Jdc Kdc] =  findsaddle(Vdc,X,Y,Z,3,Zval);
%warn('DC',Vdc,Idc,Jdc,Kdc);
%plotpot(Vdc,Idc,Jdc,Kdc,grid,0,'DC potential (no stray field included)','V_{dc} (Volt)');

if ~isempty(params.E),      % check if the initial guess for E is ok
    EE = params.E;
    Vdc = CalcVDC(data,scale*VELDC,EE(1),EE(2),EE(3),NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset); 
    % DC parameters                                                            
    [Idum Jdum Kdum] =  findsaddle(Vdc,X,Y,Z,3,Zval);
%plotpot(Vdc,Idum,Jdum,Kdum,grid,dcplot,'DC potential (stray field included)','V_{dc} (Volt)',outpath,newfilename);
    if (warn('DC',Vdc,Idum,Jdum,Kdum))&&(~strcmp(func,'N')), 
        params.E = []; 
        isempty(params.E)
    end
end

% need to restore this for d_e to run properly
Vdc = CalcVDC(data,scale*VELDC,0,0,0,NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset);

if strcmp(func,'E'),                                                       
    % this option means find stray field
    %E0 = 1e-3*[-470; -750; 24];
    while 1
    if isempty(params.E)
        while 1
            st = input('Give an initial guess for stray field (in V/m).\n','s');
            E0 = sscanf(st,'%f',inf)'/1e3;
            dist0 = d_e(E0);
            %Vdum = VDC1(scale*(W-Hor),scale*(N+Hor),scale*(Cnt+Ver),E0(1),E0(2),E0(3));
            Vdum = CalcVDC(data,scale*VELDC,E0(1),E0(2),E0(3),NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset);
            [Idum Jdum Kdum] =  findsaddle(Vdum,X,Y,Z,3,Zval);
            warn('DC',Vdum,Idum,Jdum,Kdum);
            plotpot(Vdum,Irf,Jrf,Krf,data.grid,2,'Initial guess for DC potential','V_{dc} (Volt)',outpath,newfilename);
            st = input('Happy (y/n)?\n','s');
            if strcmp(st,'y'), break; end
        end
    else
        E0 = params.E;
        dist0 = d_e(E0);
        %Vdum = VDC1(scale*(W-Hor),scale*(N+Hor),scale*(Cnt+Ver),E0(1),E0(2),E0(3));
        Vdum = CalcVDC(data,scale*VELDC,E0(1),E0(2),E0(3),NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset);
        [Idum Jdum Kdum] =  findsaddle(Vdum,X,Y,Z,3,Zval);
        warn('DC',Vdum,Idum,Jdum,Kdum);
        plotpot(Vdum,Idum,Jdum,Kdum,data.grid,0,'Initial guess for DC potential','V_{dc} (Volt)',outpath,newfilename);
    end
    fprintf('Initial guess for stray field: ( %G, %G, %G) V/m.\n',1e3*E0(1),1e3*E0(2),1e3*E0(3));
    fprintf('Miscompensation in the presence of this field: %G micron.\n\n',1e3*dist0);
    fprintf('Optimizing stray field value...\n')
    hnd = @d_e;
    E = fminsearch(hnd,E0,optimset('TolFun',(X(2)-X(1))/200));
    dist = d_e(E);
    fprintf('Stray field is ( %G, %G, %G) V/m.\n',1e3*E(1),1e3*E(2),1e3*E(3));
    fprintf('With this field the compensation is optimized to %G micron.\n\n',1e3*dist);
    if dist>5e-3, 
        params.E = [];
        fprintf('Miscompensation larger than 5 micron. Repeating.\n');
    else
        break;
    end
    end
elseif strcmp(func,'C'),                                                   
    % this option means find compensation voltages
    if isempty(params.E),
        st = input('What is the stray field you have (in V/m)?\n','s');
        E = sscanf(st,'%f',inf)'/1e3;
    else
        E = params.E;
    end
    while 1
        st = input('Give an initial guess for compensation values (Vert,Hor).\n','s');
        guess = sscanf(st,'%f',inf);
        VC0 = guess(1); HC0 = guess(2);
        %Vdum = VDC1(scale*(W-HC0),scale*(N+HC0),scale*(Cnt+VC0),E(1),E(2),E(3));
        Vdum = CalcVDC(data,scale*VELDC,E(1),E(2),E(3),NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset);
        [Idum Jdum Kdum] =  findsaddle(Vdum,X,Y,Z,3,Zval);
        warn('DC',Vdum,Idum,Jdum,Kdum);
        plotpot(Vdum,Idum,Jdum,Kdum,grid,2,'Initial guess for DC potential','V_{dc} (Volt)',outpath,newfilename);
        st = input('Happy (y/n)?\n','s');
        if strcmp(st,'y'), break; end
    end
    C0 = [HC0 VC0];
    %C0 = [2 0.8]; % for running through all the ulm indices
    dist0 = d_c(C0);       
    fprintf('Initial guess for compensation parameters: H = %G, V = %G V.\n',HC0,VC0);
    fprintf('Miscompensation in the presence of these parameters: %G micron.\n\n',1e3*dist0);
    fprintf('Optimizing compensation values...\n')
    hnd = @d_c;
    CF = fminsearch(hnd,C0,optimset('TolFun',(X(2)-X(1))/200));
    Hor = CF(1); 
    Ver = CF(2);
    dist = d_c(CF);
    fprintf('Compensation parameters: H = %G, V = %G V.\n',Hor,Ver);
    fprintf('With these parameters compensation optimized to %G micron.\n\n',1e3*dist);
elseif strcmp(func,'N')
    % this option means do not optimize anything, and just analyze the trap
    fprintf('Running ppt2 in plain analysis mode (no optimizations).\n');
    E = params.E;
    dist = NaN;
    dist = d_e(E);
    fprintf('Stray field is ( %G, %G, %G) V/m.\n',1e3*E(1),1e3*E(2),1e3*E(3));
    fprintf('With this field the compensation is optimized to %G micron.\n\n',1e3*dist);
else
    fprintf('\nInvalid ''func'' input option. Quiting.\n');
    return;
end


Vdc = CalcVDC(data,scale*VELDC,E(1),E(2),E(3),NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset);


[XRF YRF ZRF] = exactsaddle(data.EL_RF,X,Y,Z,2,Zval);                                  % find secular frequencies etc.
[XDC YDC ZDC] = exactsaddle(Vdc,X,Y,Z,3,Zval);
fprintf('RF saddle: (%f %f %f)\nDC saddle (%f %f %f).\n',XRF,YRF,ZRF,XDC,YDC,ZDC);

plotpot(Vdc,Irf,Jrf,Krf,data.grid,dcplot,'Compensated DC potential','V_{dc} (V)',outpath,newfilename);
[IDC JDC KDC] = findsaddle(Vdc,X,Y,Z,3,Zval);

[fx fy fz theta Depth rx ry rz xe ye ze superU] = pfit(E(1),E(2),E(3));

Qrf = spherharmxp(Vrf,XRF,YRF,ZRF,7,X,Y,Z);                                         % find multipole coefficients

if sqrt((XRF-XDC)^2+(YRF-YDC)^2+(ZRF-ZDC)^2)>0.008, 
    Qdc = spherharmxp(Vdc,XRF,YRF,ZRF,7,X,Y,Z); 
else
    Qdc = spherharmxp(Vdc,XDC,YDC,ZDC,7,X,Y,Z);
end    

if qcheck
    out.QualityRF = spherharmq(Vrf,Qrf,XRF,YRF,ZRF,7,X,Y,Z,'RF')
    out.QualityDC = spherharmq(Vdc,Qdc,XDC,YDC,ZDC,7,X,Y,Z,'DC')
end

Arf = 2*sqrt( (3*Qrf(8))^2+(3*Qrf(9))^2 );
Thetarf = 45*(sign(Qrf(9)))-90*atan((3*Qrf(8))/(3*Qrf(9)))/pi
Adc = 2*sqrt( (3*Qdc(8))^2+(3*Qdc(9))^2 );
Thetadc = 45*(sign(Qdc(9)))-90*atan((3*Qdc(8))/(3*Qdc(9)))/pi

%out.W = params.scale*(params.W-params.Hor);
%out.N = params.scale*(params.N+params.Hor);
%out.Center = params.scale*(params.Center+params.Ver);
%out.Hor = 0*params.Hor;
%out.Ver = 0;
%out.scale = 1;
out.E = E;
out.miscompensation = dist;
out.ionpos = [XRF YRF ZDC];
out.ionposIndex = [Irf Jrf Krf];
out.f = [fx fy fz];
out.theta = theta;
out.trapdepth = Depth/e;
out.escapepos = [xe ye ze];
out.Quadrf = 2*[Qrf(8)*3 Qrf(5)/2 Qrf(9)*6 -Qrf(7)*3 -Qrf(6)*3];
out.Quaddc = 2*[Qdc(8)*3 Qdc(5)/2 Qdc(9)*6 -Qdc(7)*3 -Qdc(6)*3];
out.Arf = Arf;
out.Thetarf = Thetarf;
out.Adc = Adc;
out.Thetadc = Thetadc;
T = [2 -2 0 0 0;...
    -2 -2 0 0 0;...
     0  4 0 0 0; ...
     0  0 1 0 0; ...
     0  0 0 1 0; ...
     0  0 0 0 1];
out.q = (1/V0)*T*out.Quadrf';
out.alpha = (2/V0)*T*out.Quaddc';
out.Error = [X(IDC)-XDC Y(JDC)-YDC Z(KDC)-ZDC]; 
if strcmp(func,'C'), 
    out.note = 'The fields useHor and useVer are the compensation parameters that ppt2 reached. They are included in W, N, and Center'
    out.useHor = Hor;
    out.useVer = Ver;
end
out.superU = superU;

%%%%%%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function f = d_e(Ei)
    % find the miscompensation distance, d_e, for the rf and dc potential 
    % given in the main program, in the presence of stray field Ei
        dm = Ei;
        E1 = dm(1); E2 = dm(2); E3 = dm(3);
        Vl = Vdc-E1*x-E2*y-E3*z;
        [xrf yrf zrf] = exactsaddle(data.EL_RF,X,Y,Z,2,Zval);
        [xdc ydc zdc] = exactsaddle(Vl,X,Y,Z,3,Zval); 
        f = sqrt((xrf-xdc)^2+(yrf-ydc)^2+(zrf-zdc)^2);
        %some debugging residues
        %if f<0.013, 
        %    fprintf(' %f %f %f %f %f %f %f %f %f %f\n',wc(1),wc(2),wc(3),wc(4),wc(5),wc(6),wc(7),wc(8),wc(9),wc(10));
        %    fprintf(' %f %f %f %f %f %f %f %f %f %f\n',nc(1),nc(2),nc(3),nc(4),nc(5),nc(6),nc(7),nc(8),nc(9),nc(10));
        %fprintf('%f %f %f\n',xrf-xdc,yrf-ydc,zdc); end
    end

    function f = d_c(C)
    % find the miscompentaion distance,  d_c, for the rf and dc potential 
    % given in the main program, in the presence of stray field 
    % E=[E(1) E(2) E(3)] defined in the main function, and the
    % compensation parameters C(1), C(2)
        dm = C;
        hc = dm(1); vc = dm(2); %ac = C(3);
        ex = E(1); ey = E(2); ez = E(3);
        wc = scale*(W-hc); nc = scale*(N+hc);  v = scale*(Cnt+vc);  
        Vl = VDC1(wc,nc,v,ex,ey,ez);
        [xdc ydc zdc] = exactsaddle(Vl,X,Y,Z,3,Zval); 
        [xrf yrf zrf] = exactsaddle(data.EL_RF,X,Y,Z,2,Zval);
        f = sqrt((xrf-xdc)^2+(yrf-ydc)^2+(zdc-zdc)^2);
        %some debugging residues
        %if f<1.0225e-5, 
        %    fprintf(' %f %f %f %f %f %f %f %f %f %f\n',wc(1),wc(2),wc(3),wc(4),wc(5),wc(6),wc(7),wc(8),wc(9),wc(10));
        %    fprintf(' %f %f %f %f %f %f %f %f %f %f\n',nc(1),nc(2),nc(3),nc(4),nc(5),nc(6),nc(7),nc(8),nc(9),nc(10));
        %    fprintf('%f\n %f %f %f\n',vc,ex,ey,ez); end
        %if f<0.07, 
        %fprintf('%f %f %f\n',xrf-xdc,yrf-ydc,zdc); end
    end


    function [fx fy fz theta Depth Xdc Ydc Zdc Xe Ye Ze superU] = pfit(e1,e2,e3)
    % find the secular frequencies, tilt angle, and position of the dc 
    % saddle point for given combined input parameters. The stray field E 
    % has been defined in the body of the main function
        
        % find dc potential
        %Vl = VDC1(w,n,cnt,ex,ey,ez);
        Vl = CalcVDC(data,scale*VELDC,e1,e2,e3,NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset);
        
        [Xdc Ydc Zdc] = exactsaddle(Vl,X,Y,Z,3,Zval);
        
        % find pseudopotential
        % 
        % Gebhard, Oct 2010:
        % changed back to calculating field numerically in ppt2 instead directly
        % with bemsolver. this is because the slow bemsolver (new version)
        % does not output EX, EY, EZ.
        %
        
        %Vrf = RFampl*data.Vrf;
        
        % grid(4) = dx, grid(5) = dy, grid(6)= dz
        [Ex,Ey,Ez] = gradient(Vrf,1e-3*grid(4),1e-3*grid(5),1e-3*grid(6));
        
        %[Ex,Ey,Ez] = gradient(Vrf,1e-3*grid(4),1e-3*grid(5),1e-3*grid(6));
        
        Esq1 = Ex.^2 + Ey.^2 + Ez.^2;
        %Esq = (RFampl*1e3*data.Erf).^2;
        %PseudoPhi = e^2*Esq/(4*Mn*mp*Omega^2); 
        
        PseudoPhi = e^2*Esq1/(4*Mn*mp*Omega^2);
      
        plotpot(PseudoPhi/e,Irf,Jrf,Krf,data.grid,plotPseudopotential,'Pseudopotential','U_{ps} (eV)',outpath,newfilename);
        
        % total trap potential
        U = PseudoPhi+e*Vl;
        
        superU = U;
        
        plotpot(U/e,Irf,Jrf,Krf,data.grid,plotTrappotential,'TrapPotential','U_{sec} (eV)',outpath,newfilename);
        %[I J K] = findsaddle(U/max(max(max(U))),X,Y,Z,3,Zval);

        
      % plotpot(data.EL_RF,Irf,Jrf,Krf,data.grid,plotTrappotential,'RF potential','(eV)',outpath,'rf');

        
        
        
        % trap frequencies and tilt in radial directions
        Uxy = U(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf); 
        MU = max(max(Uxy)); dL = (x(Irf+3,Jrf,Krf)-x(Irf,Jrf,Krf));
        Uxy = Uxy/MU;
        xr = (x(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-x(Irf,Jrf,Krf))/dL ;
        yr = (y(Irf-3:Irf+3,Jrf-3:Jrf+3,Krf)-y(Irf,Jrf,Krf))/dL;
        [C1 C2 theta] = p2d(Uxy,xr,yr);                        % look at TestMyFunctions for testiing
        
        XRF;
        x(Irf,Jrf,Krf);
        
        YRF;
        [y(Irf,Jrf-1,Krf) y(Irf,Jrf,Krf) y(Irf,Jrf+1,Krf)];
        
        
        % notes (gebhard):
        % dL = ??
        % C1, C2 coefficients of potential curve, i.e. if U = 1/2 k x^2,
        % then C1 = 1/2 k
        
        fx = (1e3/dL)*sqrt(2*C1*MU/(Mn*mp))/(2*pi); 
        fy = (1e3/dL)*sqrt(2*C2*MU/(Mn*mp))/(2*pi);
        
        % trap frequency in axial direction
        MU = 1;
        Uz = projection(U,Irf,Jrf,3)/MU;
        l1 = max([Krf-6 1]); l2 = min([Krf+6 size(Z,1)]);
        p = polyfit((Z(l1:l2)-Z(Krf))/dL,Uz(l1:l2),6); ft = polyval(p,(Z-Z(Krf))/dL,6); 
        %((Z(l1:l2)-Z(Krf))/dL)'
        %Uz(l1:l2)'
        %plot(Z,MU*Uz); hold on; plot(Z(l1:l2),MU*ft(l1:l2),'r'); 
        %title('Potential in axial direction');
        %xlabel('axial direction (mm)'); ylabel('trap potential (J)');hold off; pause(0.1);
        fz= (1e3/dL)*sqrt(2*p(5)*MU/(Mn*mp))/(2*pi);
        [Depth Xe Ye Ze] = trapdepth(U,X,Y,Z,Irf,Jrf,Krf);                
    end

    function out = normalize(in)
    % keep only the first 4 significant digits of the increment in vector
    % "in"
        dr = (in(size(in,1))-in(1))/(size(in,1)-1);
        p = 0; cnt = 0;
        while (cnt == 0)
            dr = 10*dr;
            cnt = fix(dr);
            p = p+1;
        end
        out = roundn(in,-p-4);       
    end
    
    function f= warn(str,Vi,I,J,K)
    % print a warning if the saddle point is out of bounds   
        f = false;
        if (I == 1)||(I == size(Vi,1)),
                fprintf('%s saddle out of bounds (I=%i).\n',str,I);
                f = true;
        end
        if (J == 1)||(J == size(Vi,2)),
                fprintf('%s saddle out of bounds (J=%i).\n',str,J);
                f = true;
        end
        if (K == 1)||(K == size(Vi,3)),
                fprintf('%s saddle out of bounds (K=%i).\n',str,K);
                f = true;
        end  
    end

end        
