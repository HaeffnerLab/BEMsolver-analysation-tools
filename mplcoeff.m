function M = mplcoeff(data,position,L)
% function M = mplcoeff(data,position,L)
%       ( multipoles    electrodes ->       )
% M =   (     V                             )
%       (                                   )
% Calculate the matrix of multipole coefficients for all of the dc
% electrodes. The electrodes E1, ..., E21 -> W1, ..., W10, N1, ..., N10, CNT
% Multipole coefficients only up to order 8 are kept, but the
% coefficients can be calculated up to order 11.
% The order of multipole coefficients is:
% 1/r0^2* [ (x^2-y^2)/2 (2z^2-x^2-y^2)/2 xy/2 yz/2 xz/2 ], where r0 is 1 mm
% (unless rescaling is applied)
% data is the cpo simulation data structure. 
% the calling function is tikhonov.m
% position is the axial position where the ion sits.
% L is the order of expansion
% Nikos June 2009

X = normalize(data.X); Y = normalize(data.Y); Z = normalize(data.Z);
[y x z] = meshgrid(Y,X,Z);
[Xrf Yrf Zrf] = exactsaddle(data.Vrf,X,Y,Z,2,position);
[Irf Jrf Krf] = findsaddle(data.Vrf,X,Y,Z,2,position);
plotpot(data.Vrf,Irf,Jrf,Krf,data.grid,2,'RF potential','Vrf (Volt)');
E = [0 0 0];
ord = [L L L L L L L L L L L L L L L L L L L L L];
%ord =  [5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5];
c = [ 1  0  0  0  0  0  0  0  0; ...
      0  0 -1  0  0  0  0  0  0; ...
      0  0  0 -1  0  0  0  0  0; ...
      0  1  0  0  0  0  0  0  0; ...
      0  0  0  0  0  0  0  6  0; ...
      0  0  0  0  1  0  0  0  0; ...
      0  0  0  0  0  0  0  0 12; ...
      0  0  0  0  0  0 -6  0  0; ...
      0  0  0  0  0 -6  0  0  0];
for el = 1:21
    W = [0 0 0 0 0 0 0 0 0 0]; N = W; Center = 0;
    if ~fix((el-1)/10), 
        W(mod(el,11)) = 1;
    elseif fix(el/21), 
        Center = 1; 
    else
        N(1+mod(el,11)) = 1; 
    end
    Vdc = VDC1(W,N,Center,E(1),E(2),E(3));
    plotpot(Vdc,Irf,Jrf,Krf,data.grid,0,sprintf('El. %i DC Potential',el),'V (Volt)');
        a = W
        a = N
        a = Center
    %Q = spherharmxp(Vdc(Irf-12:Irf+12,Jrf-12:Jrf+12,Krf-12:Krf+12),Xrf,Yrf,Zrf,ord(el),X(Irf-12:Irf+12),Y(Jrf-12:Jrf+12),Z(Krf-12:Krf+12));                        % this is a column [Q1 Q2 ...]'
    Q = spherharmxp(Vdc,Xrf,Yrf,Zrf,ord(el),X,Y,Z);                        % this is a column [Q1 Q2 ...]'
    M1(:,el) = Q(1:(L+1)^2);
    %M1(:,el) = Q(1:36);
end
size(M1);
M = vertcat(c*M1(1:9,:),M1(10:(L+1)^2,:)); 
%M = vertcat(c*M1(1:9,:),M1(10:36,:)); 

%%%%%%%%%%%%%%%%%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

    function Vout = VDC1(W,N,Center,EX,EY,EZ)
    % dc potential due to the cct voltage parameters and stray field EX,EY,EZ    
        Vout = W(1)*data.W1+N(1)*data.N1+W(2)*data.W2+N(2)*data.N2+W(3)*data.W3+N(3)*data.N3 ...
             +W(4)*data.W4+N(4)*data.N4+W(5)*data.W5+N(5)*data.N5+W(6)*data.W6+N(6)*data.N6 ...
             +W(7)*data.W7+N(7)*data.N7+W(8)*data.W8+N(8)*data.N8+W(9)*data.W9+N(9)*data.N9 ...
             +W(10)*data.W10+N(10)*data.N10+Center*data.Vc-EX*x-EY*y-EZ*z;
    end

end