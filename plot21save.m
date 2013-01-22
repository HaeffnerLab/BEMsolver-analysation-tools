function plot21(V,outpath,newfilename,pos)
% mesh the values of the dc voltages corresponding to the 21 dc electrodes
% of a planar trap, in a geometrically "correct" way
% V is a vector of 21 elements
% Nikos, July 2009

if max(size(V)-[21 1]),
    fprintf('plot 21: Wrong input value! Quitting\n');
    return;
end

A = zeros(100,90);
for i=1:10
    A(10*(i-1)+1:10*i,11:20) = V(i);
end
A(:,41:50) = V(21)*ones(100,10);
for i=1:10
    A(10*(i-1)+1:10*i,71:80) = V(i+10);
end
h=figure;
mesh(A);
axis([0 100 0 100 -20 20]);
print(h,'-djpeg',[sprintf('%s%s/jpeg/%d',outpath,newfilename,pos) 'Voltages3D.jpeg'])
%shading interp;