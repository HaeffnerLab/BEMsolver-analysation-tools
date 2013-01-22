function datout = trapknobs(data,position,plotOpt,superC,reg)
% function datout = trapknobs(data,position)
% return a 8x21 field datout.K with the linear combimations of trap 
% electrode voltages that give 1 V/mm, or 1 V/mm^2 of the multipole number 
% i.
% The order of multipole coefficients is:
% 1/r0 ^ [ x y z ] and 
% 1/r0^2* [ (x^2-y^2)/2 (2z^2-x^2-y^2)/2 xy/2 yz/2 xz/2 ], where r0 is 1 mm
% (unless rescaling is applied)
%

% data:     input data struct
% position: position of the ion inside the trap
% plotOpt:  show plots, 1 = yes
% superC:   if running through anlz_mlp_transport.m, this will allow for
%           forwarding previous voltages to this file (in the for-loop);
%           if you not using transport calculations, set it to e.g. 0
% reg:      sets whether regularization is done or not. 1 = yes

fprintf('executing trapknobs\n');
datout = data;
SIunits;
M = data.M;
Mt = vertcat(M(2:9,:));%,eye(21));
C = [];
c = [ 0  1  0  0  0  0  0  0; ...
      0  0  1  0  0  0  0  0; ...
     -1  0  0  0  0  0  0  0; ...
      0  0  0  0  0  0  6  0; ...
      0  0  0  1  0  0  0  0; ...
      0  0  0  0  0  0  0 12; ...
      0  0  0  0  0 -6  0  0; ...
      0  0  0  0 -6  0  0  0];
  
if isfield(data, 'allmultipoles') && data.allmultipoles == 0
    Mt(4,:) = 0;
    Mt(6,:) = 0;
    Mt(7,:) = 0;
    Mt(8,:) = 0;
end

%pause;
  
currentPos=0;
for ii=1:8
        
        Mf = zeros(8,1);
        Mf(ii) = 1;
        P = Mt\Mf;
        Mout = Mt*P; err = Mf-Mout;

        if plotOpt
            plot(err); title('Error of the fit elements');
            pause;
            plot21(P);
            pause;
        end

        %imagesc21(P,'');

        C(ii,:) = P';
        
end

K = null(Mt);


if reg
    for ii = 1:8
        Cv = C(ii,:)';
        lambda = K\Cv;
        C(ii,:) = C(ii,:)-(K*lambda)';
    end
end


% if 0
%     [Xrf Yrf Zrf] = exactsaddle(data.Vrf,data.X,data.Y,data.Z,2,position);
%     W = [0 0 0 0 0 0 0 0 0 0];
%     N = W; 
%     Center = 0;
%     for ii=1:8
%         W = W+C(ii,1:10);
%         N = N+C(ii,11:20);
%         Center = Center+C(ii,21);
%         Vii = W(1)*data.W1+N(1)*data.N1+W(2)*data.W2+N(2)*data.N2+W(3)*data.W3+N(3)*data.N3 ...
%              +W(4)*data.W4+N(4)*data.N4+W(5)*data.W5+N(5)*data.N5+W(6)*data.W6+N(6)*data.N6 ...
%              +W(7)*data.W7+N(7)*data.N7+W(8)*data.W8+N(8)*data.N8+W(9)*data.W9+N(9)*data.N9 ...
%              +W(10)*data.W10+N(10)*data.N10+Center*data.Vc;
%         Qii = spherharmxp(Vii,Xrf,Yrf,Zrf,11,data.X,data.Y,data.Z);
%         fprintf('C %i: %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f %7.6f\n',ii,c*Qii(2:9));
%     end
% end

datout.K = K;
datout.C = C';


