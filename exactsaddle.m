function [Xs Ys Zs] = exactsaddle(V,X,Y,Z,dim,Z0)
% [Xs Ys Zs] = exactsaddle(V,X,Y,Z,dim,Z0)
% This version finds the approximate saddle point using the
% pseudopotential, does a multipole expansion around it, and finds the 
% exact saddle point by maximizing the quadrupole terms
% Data are stored in matrix V. The axes are: 
% x = radial horizontal -> I
% y = radial vertical   -> J
% z = axial             -> K
% X,Y,Z are vectors defining the grid in the three directions.
% dim is the dimensionality 2, or 3
% Z0 is the coordinate where a saddle point will be sought if dim==2
% Nikos Daniilidis 9-1-09


Xs = 0; Ys = 0; Zs = 0; 

if dim==3,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize 
    %limits = [min(X) max(X) min(Y) max(Y) min(Z) max(Z)];
    %gridsize=[X(2)-X(1) Y(2)-Y(1) Z(2)-Z(1)];
    grid = [X(1) Y(1) Z(1) X(2)-X(1) Y(2)-Y(1) Z(2)-Z(1)];
    [I J K] = findsaddle(V,X,Y,Z,3,Z0);
    if (I<3)||(I>size(V,1)-2), 
        fprintf('exactsaddle.m: Saddle point out of bounds in radial direction.\n');
        return; 
    end
    if (J<3)||(J>size(V,2)-2), 
        fprintf('exactsaddle.m: Saddle point out of bounds in vertical direction.\n');
        return; 
    end
    if (K<3)||(K>size(V,3)-2), 
        fprintf('exactsaddle.m: Saddle point out of bounds in axial direction.\n');
        return; 
    end
  
    A = V(I-2:I+2,J-2:J+2,K-2:K+2); 
    xo = X(I); yo = Y(J); zo = Z(K);
    %limits = [X(I-2) X(I+2) Y(J-2) Y(J+2) Z(K-2) Z(K+2)];
    %gridsize=[(limits(2)-limits(1))/4 (limits(4)-limits(3))/4
    %(limits(6)-limits(5))/4];
    Xn = double(X(I-2:I+2)); Yn = double(Y(J-2:J+2)); Zn = double(Z(K-2:K+2));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Minimize
    
    h = @D;
    r = fminsearch(h,[xo, yo, zo]);

    Xs = r(1); 
    Ys = r(2); 
    Zs = r(3); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dim==2,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize
    if numel(size(V))==3,
        K = find(Z<Z0,1,'last');
        if K==size(Z,2), return; end;
        v1 = V(:,:,K);
        v2 = V(:,:,K+1);
        V2 = v1+(v2-v1)*(Z0-Z(K))/(Z(K+1)-Z(K));
    %limits = [min(X) max(X) min(Y) max(Y)];
    %gridsize=[X(2)-X(1) Y(2)-Y(1)];
    end
    [I J] = findsaddle(V2,X,Y,Z,2,Z0);
    if (I<3)||(I>size(V,1)-2), 
        fprintf('exactsaddle.m: Saddle point out of bounds in radial direction.\n');
        return; 
    end
    if (J<3)||(J>size(V,2)-2), 
        fprintf('exactsaddle.m: Saddle point out of bounds in vertical direction.\n');
        return; 
    end
    for k=1:9,
        A(:,:,k) = V2(I-4:I+4,J-4:J+4);
    end
    xo = X(I); yo = Y(J); z0 = Z0;
    %limits = [X(I-2) X(I+2) Y(J-2) Y(J+2) Z(K-2) Z(K+2)];
    %gridsize=[(limits(2)-limits(1))/4 (limits(4)-limits(3))/4 (limits(6)-limits(5))/4];
    Xn = X(I-4:I+4); Yn = Y(J-4:J+4); Zn = Z(K-4:K+4);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Minimize
    h = @D2;
    r = fminsearch(h,[xo, yo]);

    Xs = r(1); 
    Ys = r(2); 
    Zs = Z0; 
end


%%% Auxiliary functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function f = D(r)
    % find the weight of high order multipole terms compared to the weight 
    % of second order multipole terms in matrix V, when the center of the
    % multipoles is at x0,y0,z0
        x0 = r(1); y0 = r(2); z0 = r(3);
        %c = spherharmxp(A,x0,y0,z0,3,limits,gridsize);
        c = spherharmxp(A,x0,y0,z0,3,Xn,Yn,Zn);
        s = c.^2; 
        %f = (sum(s(2:4))+sum(s(10:16)))/sum(s(5:9));
        f = (sum(s(2:4)))/sum(s(5:9));
    end

    function f = D2(r)
    % find the weight of high order multipole terms compared to the weight 
    % of second order multipole terms in matrix V, when the center of the
    % multipoles is at x0,y0,z0
        x0 = r(1); y0 = r(2);
        %c = spherharmxp(A,x0,y0,z0,3,limits,gridsize);
        c = spherharmxp(A,x0,y0,z0,4,Xn,Yn,Zn);
        s = c.^2; 
        %f = (sum(s(2:4))+sum(s(10:16)))/sum(s(5:9));
        f = (sum(s(2:4)))/sum(s(5:9));
    end
    
end


