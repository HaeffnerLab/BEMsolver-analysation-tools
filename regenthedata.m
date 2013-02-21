function datout = regenthedata_Dtrap(data,Xcorrection,Ycorrection,position,L,truncVoltages,RF_offset)
% function datout = regenthedata(data,Xcorrection,Ycorrection,position,L)
% Regenerate the potential data for all electrodes using multipole
% expansion to order L. Also return a field datout.M, which contains the
% multipole coefficients for all electrodes.
% The electrodes E1, ..., E21 -> W1, ..., W10, N1, ..., N10, CNT
%       ( multipoles    electrodes ->       )
% M =   (     V                             )
%       (                                   )
% Multipole coefficients only up to order 8 are kept, but the
% coefficients of CNT are calculated up to order L.
% data is the cpo simulation data structure. 
% Xcorrection and Ycorrection are correction offsets from the RF saddle
% point, in case that was wrong
% position is the axial position where the ion sits.
% L is the order of expansion
% Nikos June 2009

% the correction Xrf,Yrf are parameters allowing one to offset the RF
% saddle point, for example to fix a wrong RF simulation
fprintf('Correction of XRF: %f mm.\n',Xcorrection);
fprintf('Correction of YRF: %f mm.\n',Ycorrection);

datout = data;
X = normalize(data.X); Y = normalize(data.Y); Z = normalize(data.Z);
[y x z] = meshgrid(Y,X,Z);

ord = zeros(1,data.NUM_DC+data.NUM_CENTER);
ord(:)=L;

Irf = 15; Jrf = 15; Krf = 15;
Xrf = X(Irf); Yrf=Y(Jrf); Zrf=Z(Krf);
Qrf = spherharmxp(data.EL_RF,Xrf,Yrf,Zrf,9,X,Y,Z);  
eval(sprintf('datout.%s = spherharmcmp(Qrf,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,L,X,Y,Z);','EL_RF'));


[Xrf Yrf Zrf] = exactsaddle(data.EL_RF,X,Y,Z,2,position);
[Irf Jrf Krf] = findsaddle(data.EL_RF,X,Y,Z,2,position);


Qrf = spherharmxp(data.EL_RF,Xrf,Yrf,Zrf,L,X,Y,Z);  

datout.Qrf = 2*[Qrf(8)*3 Qrf(5)/2 Qrf(9)*6 -Qrf(7)*3 -Qrf(6)*3];
datout.thetarf = 45*(sign(Qrf(9)))-90*atan((3*Qrf(8))/(3*Qrf(9)))/pi;

E = [0 0 0];

%ord = [L L L L L L L L L L L L L L L L L L L L L];



%ord =  [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
c = [ 1  0  0  0  0  0  0  0  0; ...
      0  0  1  0  0  0  0  0  0; ...
      0  0  0  1  0  0  0  0  0; ...
      0 -1  0  0  0  0  0  0  0; ...
      0  0  0  0  0  0  0  6  0; ...
      0  0  0  0  1  0  0  0  0; ...
      0  0  0  0  0  0  0  0 12; ...
      0  0  0  0  0  0 -6  0  0; ...
      0  0  0  0  0 -6  0  0  0];



if (truncVoltages == 0)       %% For a regular un-constrained trap
      for el = 1:(data.NUM_DC+data.NUM_CENTER)
        VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
        VELDC(el) = 1;

        Vdc = CalcVDC(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,2,sprintf('El. %i DC Potential',el),'V (Volt)');
%         position
%         el
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
      end

elseif (truncVoltages == 1) %% If the electrodes are connected in pairs.

    %here we have to define the electrodes that are  INCLUDED in the
    %calculation, i.e. here we take into accoung 1:8, 21:28 and 41:42.
    for el = [2 4 6 8 10 12 14 16 22 24 26 28 30 32 34 36 41:42]
        VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
        VELDC(el) = 1;
        % this connects electrodes in pairs of two (except center
        % electrodes 41,42)
        if ~(el==42||el==41) 
            VELDC(el) = 1;
            VELDC(el+1) = 1;
        end;

        Vdc = CalcVDC_Dtrap(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,2,sprintf('El. %i DC Potential',el),'V (Volt)');
%         position
%         el
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
    end   
    % resort the matrix:
    M12(:,1) = M1(:,2);
    M12(:,2) = M1(:,4);
    M12(:,3) = M1(:,6);
    M12(:,4) = M1(:,8);
    M12(:,5) = M1(:,10);
    M12(:,6) = M1(:,12);
    M12(:,7) = M1(:,14);
    M12(:,8) = M1(:,16);
    %M12(:,9) = M1(:,9);

    M12(:,9) = M1(:,22);
    M12(:,10) = M1(:,24);
    M12(:,11) = M1(:,26);
    M12(:,12) = M1(:,28);
    M12(:,13) = M1(:,30);
    M12(:,14) = M1(:,32);
    M12(:,15) = M1(:,34);
    M12(:,16) = M1(:,36);
    %M12(:,18) = M1(:,29);
    
    M12(:,17) = M1(:,41);
    M12(:,18) = M1(:,42);

    M1 = M12;
% no truncation, it will take care of the whole trap, prepare the multipole
% matrix for all electrodes.
elseif(truncVoltages == 2)
   %% For Squip first trail
     for el = [1:8 12:19 23]    % Trapping with fewer electrodes 
        VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
        VELDC(el) = 1;

        Vdc = CalcVDC_Dtrap(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,2,sprintf('El. %i DC Potential',el),'V (Volt)');
        position
        el
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
     end
% %          resort the matrix:
%     M12(:,1) = M1(:,1);
%     M12(:,2) = M1(:,2);
%     M12(:,3) = M1(:,3);
%     M12(:,4) = M1(:,4);
%     M12(:,5) = M1(:,5);
%     M12(:,6) = M1(:,6);
%     M12(:,7) = M1(:,7);
%     M12(:,8) = M1(:,8);
% %     M12(:,9) = M1(:,9);
% 
%     M12(:,9) = M1(:,12);
%     M12(:,10) = M1(:,13);
%     M12(:,11) = M1(:,14);
%     M12(:,12) = M1(:,15);
%     M12(:,13) = M1(:,16);
%     M12(:,14) = M1(:,17);
%     M12(:,15) = M1(:,18);
%     M12(:,16) = M1(:,19);
% %     M12(:,17) = M1(:,20);
% % 
%     M12(:,17) = M1(:,23);
%   
%     M1 = M12;
elseif (truncVoltages == 3)
     for el = [2:9 13:20 23]    % Trapping with fewer electrodes 
        VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
        VELDC(el) = 1;

        Vdc = CalcVDC_Dtrap(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,2,sprintf('El. %i DC Potential',el),'V (Volt)');
        position
        el
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
     end
%          resort the matrix:
    M12(:,1) = M1(:,2);
    M12(:,2) = M1(:,3);
    M12(:,3) = M1(:,4);
    M12(:,4) = M1(:,5);
    M12(:,5) = M1(:,6);
    M12(:,6) = M1(:,7);
    M12(:,7) = M1(:,8);
    M12(:,8) = M1(:,9);
%     M12(:,9) = M1(:,9);

    M12(:,9) = M1(:,13);
    M12(:,10) = M1(:,14);
    M12(:,11) = M1(:,15);
    M12(:,12) = M1(:,16);
    M12(:,13) = M1(:,17);
    M12(:,14) = M1(:,18);
    M12(:,15) = M1(:,19);
    M12(:,16) = M1(:,20);
%     M12(:,17) = M1(:,20);
% 
    M12(:,17) = M1(:,23);
  
    M1 = M12;
 
elseif (truncVoltages == 35)
     for el = [2:9 13:20 23]    % Trapping with fewer electrodes 
        VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
        VELDC(el) = 1;

        Vdc = CalcVDC_Dtrap(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,2,sprintf('El. %i DC Potential',el),'V (Volt)');
        position
        el
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
     end

elseif (truncVoltages == 5)
%% Trapping at 3rd electrode with 13 electrodes   
    for el = [1:5 12:16 23]    % Trapping with fewer electrodes 
        VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
        VELDC(el) = 1;

        Vdc = CalcVDC_Dtrap(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,2,sprintf('El. %i DC Potential',el),'V (Volt)');
        position
        el
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
     end
% %          resort the matrix:
    M12(:,1) = M1(:,1);
    M12(:,2) = M1(:,2);
    M12(:,3) = M1(:,3);
    M12(:,4) = M1(:,4);
    M12(:,5) = M1(:,5);
%     M12(:,6) = M1(:,6);
%     M12(:,7) = M1(:,7);
%     M12(:,8) = M1(:,8);
%     M12(:,9) = M1(:,9);

    M12(:,6) = M1(:,12);
    M12(:,7) = M1(:,13);
    M12(:,8) = M1(:,14);
    M12(:,9) = M1(:,15);
    M12(:,10) = M1(:,16);
%     M12(:,14) = M1(:,17);
%     M12(:,15) = M1(:,18);
%     M12(:,16) = M1(:,19);
%     M12(:,17) = M1(:,20);
% 
    M12(:,11) = M1(:,23);
  
    M1 = M12;

%trunc=6: For SQIP Palladium trap 17 electrodes| paired electrodes:1,2|9,10,11|12,13|20,21,22 | 2/15/2012
elseif (truncVoltages == 6)
        %% For sqip with some of the electrodes paired    
       for el = [1 3 4 5 6 7 8 9 12 14 15 16 17 18 19 20 23]
        VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
        VELDC(el) = 1;
        % this connects some of the electrodes pairwise
        if (el==1 || el==12) 
            VELDC(el) = 1;
            VELDC(el+1) = 1;
        end;
        if (el==9 || el==20)
            VELDC(el) = 1;
            VELDC(el+1) = 1;
            VELDC(el+2) = 1;
        end;
        disp(el)
        disp('--------------')
        disp(VELDC);

        Vdc = CalcVDC_Dtrap(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,2,sprintf('El. %i DC
        %Potential',el),'V (Volt)');
%         position
%         el
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
    end   
    % resort the matrix:
    M12(:,1) = M1(:,1);
    M12(:,2) = M1(:,3);
    M12(:,3) = M1(:,4);
    M12(:,4) = M1(:,5);
    M12(:,5) = M1(:,6);
    M12(:,6) = M1(:,7);
    M12(:,7) = M1(:,8);
    M12(:,8) = M1(:,9);
    %M12(:,9) = M1(:,9);

    M12(:,9) = M1(:,12);
    M12(:,10) = M1(:,14);
    M12(:,11) = M1(:,15);
    M12(:,12) = M1(:,16);
    M12(:,13) = M1(:,17);
    M12(:,14) = M1(:,18);
    M12(:,15) = M1(:,19);
    M12(:,16) = M1(:,20);
    %M12(:,18) = M1(:,29);
    
    M12(:,17) = M1(:,23);

    M1 = M12;
    
%trunc=7: For SQIP 200mu Gold trap 17 electrodes/ Center electrode grounded| paired electrodes:1,2|9,10,11|12,13|20,21,22 | 4/12/2012
elseif (truncVoltages == 7)
        %% For sqip with some of the electrodes paired    
       for el = [1 3 4 5 6 7 8 9 12 14 15 16 17 18 19 20]
        VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
        VELDC(el) = 1;
        % this connects some of the electrodes pairwise
        if (el==1 || el==12) 
            VELDC(el) = 1;
            VELDC(el+1) = 1;
        end;
        if (el==9 || el==20)
            VELDC(el) = 1;
            VELDC(el+1) = 1;
            VELDC(el+2) = 1;
        end;
        disp(el)
        disp('--------------')
        disp(VELDC);

        Vdc = CalcVDC_Dtrap(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        %plotpot(Vdc,Irf,Jrf,Krf,data.grid,2,sprintf('El. %i DC
        %Potential',el),'V (Volt)');
%         position
%         el
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
    end   
    % resort the matrix:
    M12(:,1) = M1(:,1);
    M12(:,2) = M1(:,3);
    M12(:,3) = M1(:,4);
    M12(:,4) = M1(:,5);
    M12(:,5) = M1(:,6);
    M12(:,6) = M1(:,7);
    M12(:,7) = M1(:,8);
    M12(:,8) = M1(:,9);
    %M12(:,9) = M1(:,9);

    M12(:,9) = M1(:,12);
    M12(:,10) = M1(:,14);
    M12(:,11) = M1(:,15);
    M12(:,12) = M1(:,16);
    M12(:,13) = M1(:,17);
    M12(:,14) = M1(:,18);
    M12(:,15) = M1(:,19);
    M12(:,16) = M1(:,20);
    %M12(:,18) = M1(:,29);
    
    %M12(:,17) = M1(:,23);

    M1 = M12;
    
    %% For the D trap MIT trapping     
elseif (truncVoltages == 4)
     for el = [6 7 8 9 17 18 19 20 23]
         VELDC = zeros(1,data.NUM_DC+data.NUM_CENTER);
         VELDC(el) = 1;
         
         Vdc = CalcVDC_Dtrap(data,VELDC,E(1),E(2),E(3),data.NUM_DC,data.NUM_CENTER,x,y,z,0,RF_offset);
        Q = spherharmxp(Vdc,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
        %eval(sprintf('datout.%s = spherharmcmp(Q,Xrf+Xcorrection,Yrf+Ycorrection,Zrf,ord(el),X,Y,Z);',str));
        M1(:,el) = Q(1:(L+1)^2);
     end
end

size(M1)
datout.M = vertcat(c*M1(1:9,:),M1(10:(L+1)^2,:)); 
size(vertcat(c*M1(1:9,:),M1(10:(L+1)^2,:)))


%%%%%%%%%%%%%%%%%%% Auxiliary functions
    
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


end