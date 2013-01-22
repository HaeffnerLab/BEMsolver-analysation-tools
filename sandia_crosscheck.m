Vsandia = load('/Users/Gebhard/Documents/Berkeley/MATLAB/data/tbStatic_0.txt');

Vsandia(2,:)=[];
Vsandia = Vsandia';
Vsandia(43:48)=[];

%Vsandia(21)=test1;
%Vsandia(22)=test2;

Vsandia2(1:20) = flipud(Vsandia(1:20));
Vsandia2(21:40) = flipud(Vsandia(22:41));
Vsandia2(41) = Vsandia(21);
Vsandia2(42) = Vsandia(42);

U = superC(:,:,5)\Vsandia2'

Vsandia2'-superC(:,:,5)*U;