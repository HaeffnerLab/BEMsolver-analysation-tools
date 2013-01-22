function printtofile2(fil,inptype,params)
% print all the results of ppt2 to file with filehandle fil. If
% fil == 1 then results are printed on the screen
% The voltages printed to the file are the actual voltages used, so they
% include the micromotion compensation value. Keep this in mind when
% calling function printtofile2
% Nikos March 2009

%fprintf(fil,'CPO simulation made on %s\n',params.datesim);
%fprintf(fil,'Voltage data taken on %s\n\n',params.datedata);

if fil ~=1,
   if strcmp(inptype,'cmb')
      V = params.cctpars(1); Df = params.cctpars(2); Cv = params.cctpars(3); A = params.cctpars(4);
      En = params.cctpars(5); H1 = params.cctpars(6); H2 = params.cctpars(7); H3 = params.cctpars(8);
      Vc = params.cctpars(9);
      fprintf(fil,'DC parameters:\nV = %G V\nDf = %G V\nCv = %G V\nA = %G V\n',V,Df,Cv,A);
      fprintf(fil,'En = %G V\nH1 = %G\nH2 = %G\nH3 = %G\nCnt = %G V\n\n',En,H1,H2,H3,Vc);
   elseif strcmp(inptype,'wsp')
      fprintf(fil,'Index of used voltages in input file: %i\n',params.index);
      fprintf(fil,'Axial ion position (this is the observed ion position): %G',params.position);
      fprintf(fil,'Compensation parameters for this index:\nHorizontal = [');
      for i=1:10
          fprintf(fil,'%G ',params.Hor(i));
      end
      fprintf(fil,']\nVertical = %G\n',params.Ver);
      fprintf(fil,'Scale = %G\n\n',params.scale);
   elseif strcmp(inptype,'ind')
       fprintf(fil,'Working with user defined input voltages');
   else
       fprintf('\nInvalid input type for printtofile.m. Quitting.\n');
       return;
   end
   fprintf(fil,'Electrode DC voltages:\nW = [')
   for i=1:10
      fprintf(fil,'%G ',params.W(i));
   end
   fprintf(fil,']\nN = [')
   for i=1:10
      fprintf(fil,'%G ',params.N(i));
   end
   fprintf(fil,']\n');
   fprintf(fil,'Cnt = %G V\n\n',params.Center);
   fprintf(fil,'RF parameters:\nVrf = %G V\nFreq = %G Hz\n\n',params.rfamplitude,params.frequency);
end

fprintf(fil,'Stray field is ( %G, %G, %G) V/m.\n',1e3*params.E(1),1e3*params.E(2),1e3*params.E(3));
fprintf(fil,'With this field the compensation is optimized to %G micron.\n',1e3*params.miscompensation);
fprintf(fil,'The ion sits at (%G,%G,%G) micron.\n',1e3*params.ionpos(1),1e3*params.ionpos(2),1e3*params.ionpos(3));
fprintf(fil,'The secular frequencies are: (%G, %G, %G) Hz.\n',params.f(1),params.f(2),params.f(3));
fprintf(fil,'The tilt angle is: %G deg.\n',params.theta);
fprintf(fil,'The trap depth is: %G eV.\n\n',params.trapdepth);

fprintf(fil,'RF potential quadrupole coefficients (in Volt/mm^2):\n');
fprintf(fil,'(x^2-y^2)/2     : %G.\n',params.Quadrf(1));
fprintf(fil,'(2z^2-x^2-y^2)/2: %G.\n',params.Quadrf(2));
fprintf(fil,'xy/2            : %G.\n',params.Quadrf(3));
fprintf(fil,'yz/2            : %G.\n',params.Quadrf(4));
fprintf(fil,'zx/2            : %G.\n',params.Quadrf(5));
fprintf(fil,'Tilted RF potential quadrupole coefficients (in Volt/mm^2):\n');
fprintf(fil,'(x''^2-y''^2)/2     : %G.\n',params.Arf);
fprintf(fil,'Tilt angle    : %G deg.\n',params.Thetarf);
fprintf(fil,'\nDC potential quadrupole coefficients (in Volt/mm^2):\n');
fprintf(fil,'(x^2-y^2)/2     : %G.\n',params.Quaddc(1));
fprintf(fil,'(2z^2-x^2-y^2)/2: %G.\n',params.Quaddc(2));
fprintf(fil,'xy/2            : %G.\n',params.Quaddc(3));
fprintf(fil,'yz/2            : %G.\n',params.Quaddc(4));
fprintf(fil,'zx/2            : %G.\n',params.Quaddc(5));
fprintf(fil,'Tilted DC potential quadrupole coefficients (in Volt/mm^2):\n');
fprintf(fil,'(x''^2-y''^2)/2     : %G.\n',params.Adc);
fprintf(fil,'Tilt angle    : %G deg.\n\n',params.Thetadc); 

A = max(params.alpha);
fprintf(fil,'\nMathieu Parameters.\n');
fprintf(fil,'alpha parameters: %G x \n',A);
fprintf(fil,'[%-5.4f %-5.4f %-5.4f]\n[%-5.4f %-5.4f %-5.4f]\n[%-5.4f %-5.4f %-5.4f]\n',params.alpha(1)/A,params.alpha(4)/A,params.alpha(6)/A, ...
                                                 params.alpha(4)/A,params.alpha(2)/A,params.alpha(5)/A, ...
                                                 params.alpha(6)/A,params.alpha(5)/A,params.alpha(3)/A);
Q = max(params.q);                                             
fprintf(fil,'q parameters: %G x\n',Q);
fprintf(fil,'[%-5.4f %-5.4f %-5.4f]\n[%-5.4f %-5.4f %-5.4f]\n[%-5.4f %-5.4f %-5.4f]\n',params.q(1)/Q,params.q(4)/Q,params.q(6)/Q, ...
                                                 params.q(4)/Q,params.q(2)/Q,params.q(5)/Q, ...
                                                 params.q(6)/Q,params.q(5)/Q,params.q(3)/Q);
%fprintf(fil,'beta parameter: %G.\n',params.betax);
%fprintf(fil,'\ny\nalpha parameter: %G.\n',params.alphay);
%fprintf(fil,'q parameter: %G.\n',params.qy);
%fprintf(fil,'beta parameter: %G.\n',params.betay);
%fprintf(fil,'\nz\nalpha parameter: %G.\n',params.alphaz);

end