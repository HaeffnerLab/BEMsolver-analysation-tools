function imagesc21(V,titl)
% function imagesc21(V)
% imagesc the values of the dc voltages corresponding to the 21 dc electrodes
% of a planar trap, in a geometrically "correct" way, and label them
% V is a vector of 21 elements
% Nikos, July 2009

if max(size(V)-[21 1]),
    fprintf('plot 21: Wrong input value! Quitting\n');
    return;
end

A = zeros(100,100);
for i=1:10
    A(10*(i-1)+1:10*i,11:30) = V(i);
end
A(:,51:60) = V(21)*ones(100,10);
for i=1:10
    A(10*(i-1)+1:10*i,71:90) = V(i+10);
end
imagesc(A); title(titl);
for i=1:10
    text(17,10*(i-1)+5,sprintf('%6.5f',V(i)));
end
text(52,50,sprintf('%6.5f',V(21)));
for i=1:10
    text(77,10*(i-1)+5,sprintf('%6.5f',V(10+i)));
end
