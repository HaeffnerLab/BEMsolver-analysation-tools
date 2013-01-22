function el = setdc(data,frequency,Ex,Ey,Ez,U1,U2,U3,U4,U5,ax,az,phi,ind,reg)
% function el = setdc(data,E,U,ax,az,phi,ind)
% set the dc voltages gor all 21 dc electrodes. 
% data is the usual data structure, E is the stray field, are the 
% five multipoles (or five 0s if no independent input is set), ax and az
% are the alpha parameters in the frame aligned with the dc quadrupole, phi 
% is the rotation of the dc quadrupole wrt the rf quadrupole in degrees, 
% ind is a boolean control that sets the control parameters to be the U's 
% (true), or the alpha parameters. 
% reg is a second boolean control determining whether the 
% output is to be regularized or not (by regularization I mean minimizing 
% the norm of el with addition of vectors belonging to the kernel of 
% data.M)
% the 'data' structure is found under ..\matlab\cpo-trap-simulations\cpo-3d\pp-data
% with name such as 'qfpdata2768_30-Jul-2009.mat' 
% Nikos, July 2009

if ind
    inp = [Ex Ey Ez U1 U2 U3 U4 U5]';
    el = data.C*inp;
else
    e=1.60217646e-19; 
    mp=1.67262158e-27;
    r0 = 1;                                                                    % length scale of multipole expansion in millimeters
    V0 = 40*mp*(2*pi*frequency)^2*(r0*1e-3)^2/e;                        % god given voltage in SI                                            % truncate & append 21x21 unity
    U2 = az*V0/8;
    %U1p = U2+ax*V0/4;
    %U1 = U1p*cos(2*pi*(phi+data.thetarf)/180);
    %U3 = 2*U1p*sin(2*pi*(phi+data.thetarf)/180);
    U1 = U2+ax*V0/4;
    U3 = 2*U1*tan(2*pi*(phi+data.thetarf)/180);
    U1p= sqrt(U1^2+U3^2/2);
    U4 = U1p*data.Qrf(4)/data.Qrf(1);
    U5 = U1p*data.Qrf(5)/data.Qrf(1);
    inp = [Ex Ey Ez U1 U2 U3 U4 U5]';
    el = data.C*inp;
end

if reg
    c = el;
    lambda = data.K\c;
    el = el-(data.K*lambda);
end
    