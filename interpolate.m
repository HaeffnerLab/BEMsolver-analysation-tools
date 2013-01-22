startInt = -0.01;
endInt = 0.1;
stepInt = 0.005;

y1 = superEl(1,:);
y2 = superEl(2,:);
y3 = superEl(3,:);

y7 = superEl(7,:);

 x = startInt:stepInt:endInt;
% 
% pp = interp1(x,y,'cubic','pp');
% y2 = ppval(pp,x);
% 
% plot(x,y,'ro'),hold on, plot(x,y2,'-'),hold off; 