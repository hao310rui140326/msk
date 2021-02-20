clear all;
close all;
clc;

FLAG_ADD_TIME   = 1 ; 
FLAG_ADD_FREQ   = 0 ;
FLAG_ADD_NOISE  = 1 ; 
flag_shift = 0  ;


shift = 0; 
EbN0 = [20]  ;%4:0.2:12;
freqOffset = 20e3;% 20e3 ; % 10e3;
ppm = 400e-6;%400e-6;%-400e-6; % residual ppm should not be too much

flag_plot = 0 ;




numFrames = 1;
samplesPerSymbol= 12;
Ns = samplesPerSymbol ;
symbolsPerFrame = 16 ;
Nf = symbolsPerFrame ;
fSample = 192e6;
%fSample = 32e6;
fSymbol = 16e6 ;
fs = fSample ;
Ts = 1.0/fs ;









    %------------------------------------------------------------------------
    % Generate and modulate data
    %------------------------------------------------------------------------
    %txBits = randi([0 1],symbolsPerFrame,1);
    err_num = zeros(length(EbN0),1) ;
    sum_err_num = zeros(length(EbN0),1) ;
    total_bits = zeros(length(EbN0),1) ;
    err_rate = zeros(length(EbN0),1) ;
    Lall = [800:100:1600]*12;   %% about 20 us*16M = 1600points
for iter = 1:100
for ll = 1:length(Lall)
    L = Lall(ll) ;
for kk = 1:1%length(ppm)
for ii = 1:length(EbN0)
for jj = 1:1 %1e6 % frame iteration
    snr = EbN0(ii) ;
    chAWGN = comm.AWGNChannel(...
    'NoiseMethod','Signal to noise ratio (Eb/No)', ...
    'EbNo',snr,...
    'SignalPower',1, ...
    'SamplesPerSymbol',samplesPerSymbol);
    hLen = 5000 ;
    %zerobits = zeros(1,1000) ;
    headBits = randi([0 1],1,3);
    tailBits = randi([0 1],1,0);
    txBits1 = [0 1 0 0 0 1 0 0 1 1 0 0 1 1 0 1] ;
    txBits = mskmod([headBits txBits1 tailBits],samplesPerSymbol) ;
    
    dLen = 32;
    txref = randi([0 1],1,dLen) ; %[ 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1] ;
    dataBits  = [ txref randi([0 1],1,144)];
    txSymData = mskmod(dataBits,samplesPerSymbol);
    
    
    txSym_t = [txBits zeros(1,L) txSymData] ;
    shift = 50 ;%randi(50) ; %% 1~50  %% ensure first point not drift too much 
    if flag_shift
        txSym = [zeros(1,shift*12) txSym_t];
    else
        txSym = txSym_t ;
    end
    %------------------------------------------------------------------------
    % Transmit through channel
    %------------------------------------------------------------------------
    if FLAG_ADD_NOISE
        rxchan = chAWGN(txSym);
    else
        rxchan = txSym ; 
    end    
    %
    % Add timing offset
    if FLAG_ADD_TIME
        Frac = 1e5; 
        timingOffset = round(ppm(kk)*Frac) ;
        rxSigTimingOff1 = resample(rxchan,Frac,Frac+timingOffset) ;
        rxSigTimingOff = rxSigTimingOff1(1:floor(length(rxSigTimingOff1)/12)*12) ;
        %     varDelay = dsp.VariableFractionalDelay;
        %     rxSigTimingOff = varDelay(rxchan,0.2);
        %pwelch(rxchan,[],[],[],192e6,'twosided');  
        if flag_plot
            pwelch(rxSigTimingOff,[],[],[],192e6,'twosided');
            stem(real(rxchan(1:192)),'b');
            hold on;
            stem(real(rxSigTimingOff(1:192)),'r') ;
        end    
    else
        rxSigTimingOff = rxchan ; 
    end    
    %
    % Add carrier frequency and phase offset
    %rxSigCFO = pfo(rxSigTimingOff);
    if FLAG_ADD_FREQ
        Len = length(rxSigTimingOff) ;
        rxSigCFO = rxSigTimingOff.*exp(-1*j*2*pi*freqOffset*[0:Len-1]*Ts) ;
    else
        rxSigCFO = rxSigTimingOff ;
    end    
    rxSig = rxSigCFO ;
    %
%     scatterplot(txSym(1:samplesPerSymbol:end));
%     scatterplot(rxSigCFO(1:samplesPerSymbol:end)) ;
    % Pass the signal through an AWGN channel
    
    %%%get the corr for ppm estimation
    %for tt=1:1
         find_loc_flag = 1 ;
         if find_loc_flag
             location =  find_location(rxSig,txBits1,16);
             lll = location(1)-1+L ;
         else
             lll = 3*12+L ;
         end    
        rxCorr = rxSig(16*12+lll+1:16*12+lll+Ns*dLen) ;  
        rxpre = rxSig(16*12+lll-23:(16*12+lll)) ;
        rxpost = rxSig(16*12+lll+1+Ns*dLen:16*12+lll+24+Ns*dLen) ;  
        
        %tx1 = txSymPr(1:192) ;
%         plot(real(tx1),'g');
%         hold on; plot(imag(tx1),'b');
%         plot(angle(tx1),'r');
%         hold on;
%         plot(angle(rxCorr),'b');
    
%     ppm_est = time_est(tx1,Ns) ;
      fprintf("current ppm is %f\n" , ppm(kk)) ;
      
      %%ideal rm should be gm(t) = 0.5*(-1)^m*(1+cos(2*pi*t/T))
      N=Ns ;
       
        % simple logic for ppm estimation
%       sum_ppm = 0 ;
%       for ii = 1:N
%           rmi = -0.5*(1+cos(2*pi*(ii-1)*(1+ppm(kk))/N)) ;
%           %fprintf("%d ideal abs %f is %f + j * %f \n" , ii-1 , abs(rmi) , real(rmi) , imag(rmi)) ;
%           sum_ppm = sum_ppm + abs(rmi)* exp(-1*j*2*pi*(ii-1)/N) ;
%       end
%       theta = angle(sum_ppm) ; 
%       t_est = -1  * theta/(2*pi) *(-2); 
%       fprintf("ppm estimation is %f, ideal %f ,mse %f \n, ",...
%             t_est , ppm(kk) , (abs(t_est-ppm(kk))/ppm(kk))^2) ;
%       break; 


        delay_est = time_est(rxCorr,Ns)*(-1) ;
        delay_ideal = ppm(kk)*(19+L/12) ;
        %scatterplot(rxCorr)
        %scatterplot(rxchan)
        fprintf("delay estimation is %f, ideal %f , residual %f ratio %f , mse %f \n",...
               delay_est , delay_ideal , abs(delay_est-delay_ideal) , (abs(delay_est)/delay_ideal) ,  (abs(delay_est-delay_ideal)/delay_ideal)^2) ;
        %keyboard;
        %%resampleing to recover original signals
        %%use lagrange interpolation
        xh = [0:Ns*dLen-1]+delay_est*Ns ;
        xx = [-2*Ns:Ns*dLen-1+2*Ns];
        yy = [rxpre rxCorr rxpost];
        yyi = zeros(1,length(xh)) ;
%         rxinterp = zeros(length(yy),1);
         for ii = 1:length(xh)
            xhi = xh(ii) ;
            ddx =[ floor(xhi)-1  floor(xhi)  floor(xhi)+1  floor(xhi)+2 ] ;
            ddy = yy(ddx+2*Ns+1) ;
            yyi(ii) =  lagrange(ddx,ddy,xhi);
         end   
        
        %% freq offset estimation 
        f_est = freq_est(yyi,Ns) ;
        fprintf("freq ideal is %f, est %f\n" ,freqOffset , f_est*fSymbol) ;
        %% residual phase estimation
        
        
        
        %% demod
        rx_demod=mskdemod(yyi,12);
        disp(rx_demod);
        rx_demod_ref=mskdemod(rxCorr,12);
                disp(rx_demod_ref);
                        disp(txref);


   % end
    %
    % Save the transmitted signal for plotting
    plot_rx = rxSig;
    break;
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	rx_demod=mskdemod(rxSig,12);
    [err_num(ii) tt] = biterr(rx_demod(1:end-1),txBits(1:end-1)) ;
    sum_err_num(ii) = sum_err_num(ii)+err_num(ii);
    total_bits(ii) = total_bits(ii) + length(txBits) ;
    
    fprintf("EbN0 is %f and FRAME NUMBER is %d\n" , snr , jj ) ;
    if( sum_err_num(ii) > 100 ) || (jj>1000) % at least 10-5
        break ;
    end  
end    
      
%     err_rate(ii) = sum_err_num(ii)/total_bits(ii) ;
%     ideal_rate(ii)=1/2*erfc(sqrt(10.^(snr/10)/2));
end
end
end
end

if flag_plot
    semilogy(EbN0,err_rate);
    hold on ;
    semilogy(EbN0,ideal_rate,'r');
end
keyboard ;

