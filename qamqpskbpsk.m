Eb = 1;
SNR_dB = 0:2:13; % in dB
SNR = 10.^(SNR_dB./10);

%number of symbols for modulation schemes
sym_BPSK = 2;
sym_QPSK = 4;
sym_QAM = 16;   
M = 4000;

%transmit-symbols for each modulation scheme
Es_BPSK = Eb*log2(sym_BPSK);
Es_QPSK = Eb*log2(sym_QPSK);
Es_QAM = Eb*log2(sym_QAM);

%transmit signal
signal_BPSK = randi([0,1],1,1000); %symbols transmitted
signal_QPSK = randi([0,3],1,2000);
signal_QAM = randi([0,15],1,1000); %16-QAM
symbol_QAM = [-3,-1,1,3];

totalerrorcount_BPSK = zeros(1,7);
totalerrorcount_QPSK = zeros(1,7);
totalerrorcount_QAM = zeros(1,7);

for MCruns = 1:15
    
    txsignal_BPSK = zeros(1,1000);
    y_BPSK = zeros(1,1000);
    txsignal_QPSK = zeros(1,2000);
    y_QPSK = zeros(1,1000);
    txsignal_QAM = zeros(1,4000);
    y_QAM = zeros(1,1000);
    
    errorcount_BPSK = zeros(1,7);
    errorcount_QPSK = zeros(1,7);
    errorcount_QAM = zeros(1,7);
    
    
    
    for n=1:1000
        
     %BPSK   
        if(signal_BPSK(n) == 0)
            txsignal_BPSK(n) = 1*exp(1i*(2*pi));
        elseif(signal_BPSK(n) == 1)
            txsignal_BPSK(n) = 1*exp(1i*(pi));
        end
        
     %QPSK
        if (signal_QPSK(n) == 0)
            txsignal_QPSK(n) = sqrt(Es_QPSK)*exp(1i*pi/4);
        elseif(signal_QPSK(n) == 1)
            txsignal_QPSK(n) = sqrt(Es_QPSK)*exp(1i*(3*pi/4));
        elseif(signal_QPSK(n) == 2)
            txsignal_QPSK(n) = sqrt(Es_QPSK)*exp(1i*(7*pi/4));
        elseif(signal_QPSK(n) == 3)
            txsignal_QPSK(n) = sqrt(Es_QPSK)*exp(1i*(5*pi/4));
        end   
                     
       
     end
     
%QAM generation

txsignal_QAM = (randsrc(1,M,symbol_QAM) + 1i*randsrc(1,M,symbol_QAM)) ;
     
     BPSK_demod = zeros(1,1000);
     QPSK_demod = zeros(1,2000);
     QAM_demod = zeros(1,1000);
             
     %received signal
    for l=1:length(SNR_dB)
     
        for n=1:1000
            %BPSK
             y_BPSK(n) = awgn(txsignal_BPSK(n),SNR_dB(l)); %received signal
             BPSK_demod(n) = pskdemod(y_BPSK(n), sym_BPSK);
            s=sign(y_BPSK);
             k1 = abs(BPSK_demod(n)-signal_BPSK(n));
             if (k1~= 0)
                 errorcount_BPSK(l) = errorcount_BPSK(l) + 1;
             end 
          
            %QPSK
            y_QPSK(n)= awgn(txsignal_QPSK(n),SNR_dB(l));
            QPSK_demod(n) = pskdemod(y_QPSK(n), sym_QPSK,pi/4,'gray');
            k2 = abs(QPSK_demod(n)-signal_QPSK(n));
            if (k2~=0)
                errorcount_QPSK(l) = errorcount_QPSK(l) + 1;
            end
            
           % QAM
             y_QAM(n) = awgn(txsignal_QAM(n),SNR_dB(l));
             QAM_demod(n) = qamdemod(y_QAM(n),sym_QAM,0,'gray');
             k3 = abs(QAM_demod(n)-txsignal_QAM(n));
            if (k3~=0)
                errorcount_QAM(l) = errorcount_QAM(l) + 1;
            end
            
            
         end
    end
         totalerrorcount_BPSK = totalerrorcount_BPSK + errorcount_BPSK;
         totalerrorcount_QPSK = totalerrorcount_QPSK + errorcount_QPSK;
         totalerrorcount_QAM = totalerrorcount_QAM + errorcount_QAM;
end


%symbol errors
 symerr_BPSK = totalerrorcount_BPSK/15000;
 symerr_QPSK = totalerrorcount_QPSK/15000;
 symerr_QAM = totalerrorcount_QAM/15000; 
 
  Pes_BPSK = zeros(1,7);
  Pes_QPSK = zeros(1,7);
  Pes_QAM = zeros(1,7);
  
%  theoretical value
  for j=1:7
  Pes_BPSK(j) = qfunc(sqrt(2*SNR(j)));
  Pes_QPSK(j) = 2*qfunc(sqrt(2*SNR(j)));
  Pes_QAM(j) = 3*qfunc(sqrt(0.8*SNR(j)))/sqrt(sym_QAM); %dmin = 0.8
  end
  
 %BPSK
 semilogy(SNR_dB,symerr_BPSK,'r',SNR_dB,Pes_BPSK);
 legend('simulated graph', 'theoretical graph');
 title('BPSK AWGN v/s Theoretical simulation');
 xlabel('SNR in dB');
 ylabel('probability of symbol error');
 
%QPSK
semilogy(SNR_dB,symerr_QPSK,'r',SNR_dB,Pes_QPSK);
legend('simulated graph', 'theoretical graph');
title('QPSK AWGN v/s Theoretical simulation');
xlabel('SNR in dB');
ylabel('probability of symbol error'); 

%QAM
scatterplot(txsignal_QAM);
title('16 QAM constellation');

semilogy(SNR_dB,Pes_QAM);% theoretical QAM graph
title('16-QAM theoretical symbol error');
xlabel('SNR in dB');
ylabel('probability of symbol error');

semilogy(SNR_dB,symerr_QAM); %simulated QAM

scatterplot(y_QAM); % constellation of noise added to QAM
 
%--------------------------------Rayleigh-----------------------------------------------------
 
RayleighSNRdb = 0:2:35;
RayleighSNR = 10.^(RayleighSNRdb./10);
Rayleigh_totalerrorcount_BPSK = zeros(1,18);
Rayleigh_totalerrorcount_QPSK = zeros(1,18);

for MCruns = 1:45
    Raytxsignal_BPSK = zeros(1,1000);
    Raytxsignal_QPSK = zeros(1,2000);
    Rayy_BPSK = zeros(1,1000);
    Rayy_QPSK = zeros(1,2000);
    
    Rayerrorcount_BPSK = zeros(1,18);
    Rayerrorcount_QPSK = zeros(1,18);
    
    for p = 1:1000
        alphac(p) = normrnd(0,0.5);
        alphas(p) = normrnd(0,0.5);
        alpha(p) = sqrt((alphac(p)).^2 + (alphas(p)).^2);
    end
    
    for m=1:1000
     
        %BPSK   
        if(signal_BPSK(n) == 0)
            Raytxsignal_BPSK(n) = 1;
        elseif(signal_BPSK(n) == 1)
            Raytxsignal_BPSK(n) = -1;
        end 
         
        %QPSK
         if (signal_QPSK(m) == 0)
             Raytxsignal_QPSK(m) = sqrt(Es_QPSK)*exp(1i*pi/4);
         elseif(signal_QPSK(m) == 1)
             Raytxsignal_QPSK(m) = sqrt(Es_QPSK)*exp(1i*(3*pi/4));
         elseif(signal_QPSK(m) == 2)
             Raytxsignal_QPSK(m) = sqrt(Es_QPSK)*exp(1i*(7*pi/4));
         elseif(signal_QPSK(m) == 3)
             Raytxsignal_QPSK(m) = sqrt(Es_QPSK)*exp(1i*(5*pi/4));
         end    
   end
     
     
    for l=1:length(RayleighSNRdb)
    
        for m=1:1000
           %BPSK
            y_BPSK(n) = awgn(alpha(n)*Raytxsignal_BPSK(n),RayleighSNRdb(l)); %received signal
            BPSK_demod(n) = pskdemod(y_BPSK(n),sym_BPSK,0,'gray');
            
            k1 = abs(BPSK_demod(n)-signal_BPSK(n));
            if (k1~= 0)
                Rayerrorcount_BPSK(l) = Rayerrorcount_BPSK(l) + 1;
            end 
            
           %QPSK
            y_QPSK(m) = awgn(alpha(m)*Raytxsignal_QPSK(m),RayleighSNRdb(l)); %received signal
            QPSK_demod(m) = pskdemod(y_QPSK(m),sym_QPSK,pi/4,'gray');
            s=sign(y_QPSK);
            k1 = abs(QPSK_demod(m)-signal_QPSK(m));
            if (k1~= 0)
                Rayerrorcount_QPSK(l) = Rayerrorcount_QPSK(l) + 1;
            end 
         end
     end
     
     Rayleigh_totalerrorcount_BPSK = Rayleigh_totalerrorcount_BPSK + Rayerrorcount_BPSK;
     Rayleigh_totalerrorcount_QPSK = Rayleigh_totalerrorcount_QPSK + Rayerrorcount_QPSK;
end

 %Rayleigh symbol error
 Raysymerr_BPSK =  Rayleigh_totalerrorcount_BPSK/45000;
Raysymerr_QPSK =  Rayleigh_totalerrorcount_QPSK/45000;

   
 semilogy(RayleighSNRdb,Raysymerr_BPSK,'b',SNR_dB,symerr_BPSK,'r');
 legend('Rayleigh graph', 'AWGN graph');
 title('BPSK Rayleigh v/s AWGN simulation');
 xlabel('SNR in dB');
 ylabel('probability of symbol error');     
   
semilogy(RayleighSNRdb,Raysymerr_QPSK,'b',SNR_dB,symerr_QPSK,'r');
legend('Rayleigh graph', 'AWGN graph');
title('QPSK Rayleigh v/s AWGN simulation');
xlabel('SNR in dB');
ylabel('probability of symbol error');     
  
   
