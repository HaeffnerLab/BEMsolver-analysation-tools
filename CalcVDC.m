function Vout = CalcVDC_Dtrap(data,VELDC,EX,EY,EZ,NUM_DC,NUM_Center,x,y,z,truncVoltages,RF_offset)

    Vout = zeros(size(data.EL_DC1));

    
    % if we want keep electrodes at the same voltages in pairs...
    % here for electrode 1 disconnected and electrodes 2-17 and 22-37
    % and 41,42 connected in pairs of two (i.e. 2,3 and 4,5 etc have same
    % voltage.
    if (truncVoltages==1)
        VELDC2 = zeros(NUM_DC+NUM_Center,1);
        for ii=1:8
            kk = ii*2;
            VELDC2(kk) = VELDC(ii);
            VELDC2(kk+1) = VELDC(ii);
        end
        % as we want to keep both electrodes at the same voltage (they are
        % shorted!!)
        %VELDC2(9) = VELDC(8);
        for ii=1:8
            kk = 20+ii*2
            VELDC2(kk) = VELDC(8+ii);
            VELDC2(kk+1) = VELDC(8+ii);
        end
        for ii=1:2
            VELDC2(40+ii) = VELDC(16+ii);
        end
        VELDC = VELDC2
    end    
    
    %% In case of only a part of the trap being used 
    if (truncVoltages==2)
        VELDC2 = zeros(NUM_DC+NUM_Center,1);
        for ii = 1:8
            VELDC2(ii) =VELDC(ii); 
            VELDC2(ii+11) =VELDC(ii+8); 
        end
        VELDC2(23) = VELDC(17);
         
        VELDC = VELDC2;
    end
    
%% Electrodes paired     
    if (truncVoltages==3)
        VELDC2 = zeros(NUM_DC+NUM_Center,1);
        for ii = 1:8
            VELDC2(ii+1) =VELDC(ii); 
            VELDC2(ii+12) =VELDC(ii+8); 
        end
        VELDC2(23) = VELDC(17);
         
        VELDC = VELDC2;
    end
%% Electorde pairing for SQIP

    if (truncVoltages==6)

        VELDC2 = zeros(NUM_DC+NUM_Center,1);
        ii = 1;
        for el = [1 3 4 5 6 7 8 9 12 14 15 16 17 18 19 20 23]
            VELDC2(el) =VELDC(ii); 
            if (el==1 ||  el==12 ) 
                VELDC2(el+1) =VELDC(ii); 
            end;
            
            if (el==9 || el==20)
                VELDC2(el+1) = VELDC(ii);
                VELDC2(el+2) = VELDC(ii);
            end;
            ii = ii+1;
        end
        VELDC = VELDC2;
    end
    
        %% In case of only a part of the trap being used 
    if (truncVoltages==5)
        VELDC2 = zeros(NUM_DC+NUM_Center,1);
        for ii = 1:5
            VELDC2(ii) =VELDC(ii); 
            VELDC2(ii+11) =VELDC(ii+5); 
        end
        VELDC2(23) = VELDC(11);
         
        VELDC = VELDC2;
    end
    
    %% No truncation
    for ii=1:NUM_DC
        Vout = Vout + VELDC(ii)*data.(['EL_DC' num2str(ii)]);
    end
    
    for ii=1:(NUM_Center)
        pos = NUM_DC+ii;
        Vout = Vout + VELDC(pos)*data.(['EL_CNT' num2str(ii)]);
        if isfield(data, 'EL_FING1')
            Vout = Vout + VELDC(pos)*data.(['EL_FING' num2str(ii)])/2;
        end
    end
    
    
    %For applying a DC offset on RF
    Vout = Vout + RF_offset*data.EL_RF;
    
    Vout = Vout-EX*x-EY*y-EZ*z;

end
