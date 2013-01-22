% clear;
% clc;
% 
% load Vset;

gap = 10;
width_DC = 600;
length_DC(1) = 600;
length_DC(2) = 200;
length_DC(3) = 900;

width_RFw = 540;
length_RFw = 15500;

width_RFn = 270;
length_RFn = 15500;

width_RFd = 180;
length_RFd = 740;

width_cen = 180;
length_cen = 14840;
ion_height = 207;

DC_sizes = [1 1 1 1 2 1 1 1 1 3 3 1 1 1 1 1 2 1 1 1 1 3 3];
DCset = Vset;
Org = [width_DC+width_RFn+width_cen/2+2*gap, 0, 4300];
I = 1;

%%%%%%%%   RF and Center Electorde Definition
    Cen1 = [width_DC+2*gap+width_RFn, 0,  length_RFd] - Org;
    Cen2 = [width_DC+2*gap+width_RFn+width_cen, 0, length_RFd] - Org;
    Cen3 = [width_DC+2*gap+width_RFn+width_cen, 0, length_RFd+length_cen] - Org;
    Cen4 = [width_DC+2*gap+width_RFn, 0, length_RFd+length_cen] - Org;
    
    RFw1 = [width_DC+3*gap+width_RFn+width_cen, 0, 0] - Org;
    RFw2 = [width_DC+3*gap+width_RFn+width_cen+width_RFw, 0, 0] - Org;
    RFw3 = [width_DC+3*gap+width_RFn+width_cen+width_RFw, 0, length_RFw] - Org;
    RFw4 = [width_DC+3*gap+width_RFn+width_cen, 0, length_RFw] - Org;
    
    RFn1 = [width_DC+2*gap, 0, 0] - Org;
    RFn2 = [width_DC+2*gap+width_RFn, 0, 0] - Org;
    RFn3 = [width_DC+2*gap+width_RFn, 0, length_RFn] - Org;
    RFn4 = [width_DC+2*gap, 0, length_RFn] - Org;
    
    RFd1 = [width_DC+2*gap+width_RFn, 0, 0] - Org;
    RFd2 = [width_DC+4*gap+width_RFn+width_RFd, 0, 0] - Org;
    RFd3 = [width_DC+4*gap+width_RFn+width_RFd, 0, length_RFd] - Org;
    RFd4 = [width_DC+2*gap+width_RFn, 0, length_RFd] - Org;

    xext = 20;         %No of array elements
    yext = 20;
    zext = 1;
    for xa = (1:xext)
        xa
        for ya = (1:yext)
            for za = (1:1)
                dx = 20.000;
                dy = 20.000;
                dz = 0.000;
                xaa = xa*dx;
                yaa = ya*dy;
                zaa = za*dz;
                RFnull = [-30, 207, 3800+6.5*gap];  %is the exact center of the RF Null for Z position = 2250
                rx = RFnull(1) - xext*dx/2;
                ry = RFnull(2) - yext*dy/2;
                rz = RFnull(3) - zext*dz/2;
%                 X = [-30+xaa, 207+yaa, 3800];  %is the exact center of the RF Null for Z position = 2250
                X = [rx+xaa, ry+yaa, rz+zaa];  %is the exact center of the RF Null for Z position = 2250
                x(xa,ya,za) = X(1) - RFnull(1);
                y(xa,ya,za) = X(2) - RFnull(2);
                z(xa,ya,za) = X(3) ;

%% RF Electric field estimation               
                Irf = -1;
                E_RF =  E_field(Irf, X, RFw1, RFw2) + E_field(Irf, X, RFw2, RFw3) + E_field(Irf, X, RFw3, RFw4) + E_field(Irf, X, RFw4, RFw1) +...
                               E_field(Irf, X, RFn1, RFn2) + E_field(Irf, X, RFn2, RFn3) + E_field(Irf, X, RFn3, RFn4) + E_field(Irf, X, RFn4, RFn1) +...
                               E_field(Irf, X, RFd1, RFd2) + E_field(Irf, X, RFd2, RFd3) + E_field(Irf, X, RFd3, RFd4) + E_field(Irf, X, RFd4, RFd1) ;
                ERFx(xa,ya) = E_RF(1);
                ERFy(xa,ya) = E_RF(2);
                ERFmag(xa,ya) = sqrt(E_RF(1)^2 + E_RF(2)^2);
            
%% DC Electric field estimation.
            E_DCtot = [0, 0, 0];
            len_offset = 4300;
                for a = 1:22
                    I = DCset(a);
                    
                    if (a>11)
                        width_add = width_DC+width_RFn+width_cen+width_RFw+4*gap;
                        elno = a-11;
                        len_offset = len_offset + length_DC(DC_sizes(a)) + gap;
                    else
                        width_add = 0;
                        elno = a;
                        len_offset = len_offset + length_DC(DC_sizes(a)) +gap;
                        len_left = len_offset;
                    end
                    %%%%%%%%a th DC electrode definition 
                    pos1 = [width_add, 0, len_offset-length_DC(DC_sizes(a))] - Org;
                    pos2 = [width_add+width_DC, 0, len_offset-length_DC(DC_sizes(a))] - Org;
                    pos3 = [width_add+width_DC, 0, len_offset] - Org;
                    pos4 = [width_add, 0, len_offset] - Org;

                    %%%%%%%%%Estimation of the eletric field for that DC electrode
                    E_DC(a,:) = E_field(I, X, pos1, pos2) + E_field(I, X, pos2, pos3)+E_field(I, X, pos3, pos4)+E_field(I, X, pos4, pos1);
                    E_DCtot = E_DCtot + E_DC(a,:);
                    
                    if (a==11)
                        len_offset = 4300;
                    end
    
                end    
                I = DCset(23);
                Efieldcen = E_field(I, X, Cen1, Cen2)+E_field(I, X, Cen2, Cen3)+E_field(I, X, Cen3, Cen4)+E_field(I, X, Cen4, Cen1); %Center electrode
                
                E_DCtot = E_DCtot + Efieldcen;
                Ex(xa,ya) = E_DCtot(1);
                Ey(xa,ya) = E_DCtot(2);
                Ez(xa,ya) = E_DCtot(3);
            end
        end
    end 
%     plot(X(1),X(3),'bo');
    
%     

figure; streamslice(x(:,1),y(1,:),ERFx',ERFy')
figure; streamslice(x(:,1),y(1,:),Ex',Ey')

figure; quiver(x,y,ERFx,ERFy)
figure; quiver(x,y,Ex,Ey)
