set(0,'DefaultAxesFontSize',14)

y=25*10^-6;
x=113*10^-6;
R=25*10^-6;
D=500*10^-6;
A=(R+D)/R
a=R/2*(A+sqrt(A^2-1)-1/(sqrt(A^2-1)+A));

Ey=-sqrt((x-a)^2+y^2)/(sqrt((x+a)^2+y^2)*((x-a)^2+y^2))*(2*y*sqrt((x-a)^2+y^2)/(sqrt((x+a)^2+y^2))-2*y*sqrt((x+a)^2+y^2)/(sqrt((x-a)^2+y^2)))/log(A+sqrt(A^2-1));
Ex=-sqrt((x-a)^2+y^2)/(sqrt((x+a)^2+y^2)*((x-a)^2+y^2))*2*((x+a)*sqrt((x-a)^2+y^2)/sqrt((x+a)^2+y^2)-(x-a)*sqrt((x+a)^2+y^2)/sqrt((x-a)^2+y^2))/log(A+sqrt(A^2-1));

dx=1/Ex*10^6;
dy=1/Ey*10^6;
dmax=1/(Ey+Ex)*10^6;

Ex2=a*log(A+sqrt(A^2-1));


datapathin = '/home/soenke/Documents/Mathlab/trap_simulations/Maxwell/wire_on_top3';

data = load([sprintf('%s',datapathin), '.reg']); 

X=data(:,1);
Y=1./abs(data(:,6)).*10^6;


figure()
plot(X,Y)
hold on
xlabel('Distance from wire [\mu m]'); ylabel('d_{eff} [\mu m]');
title('250 \mu m wire height');
axis([0,250,0,800])
hold off