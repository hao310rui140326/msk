clear all;
close all;
clc;

FLAG_ADD_TIME   = 0 ; 
FLAG_ADD_FREQ   = 1 ;
FLAG_ADD_NOISE  = 0 ; 
flag_shift = 0  ;

shift = 0; 
SNR = 20;%4:0.2:12;
%%freqOffset = -(0:0.01:0.2)*16e6 ; % 10e3;
freqOffset = -100e3 ; % 10e3;
ppm = 0 ; %0.02

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

intv0 = 500  ;
intv1 = 500  ;
intv2 = 500  ;
intv3 = 500  ;
hLen  = 5000 ;
tLen  = 5000 ;
pLen  = 300  ;

    %------------------------------------------------------------------------
    % Generate and modulate data
    %------------------------------------------------------------------------
    %txBits = randi([0 1],symbolsPerFrame,1);
    err_num = zeros(length(freqOffset),1) ;
    sum_err_num = zeros(length(freqOffset),1) ;
    total_bits = zeros(length(freqOffset),1) ;
    err_rate = zeros(length(freqOffset),1) ;
for ii = 1:length(freqOffset)
for jj = 1:1%1e6 % frame iteration
    snr = SNR ;
    chAWGN = comm.AWGNChannel(...
    'NoiseMethod','Signal to noise ratio (Eb/No)', ...
    'EbNo',snr,...
    'SignalPower',1, ...
    'SamplesPerSymbol',samplesPerSymbol);
    intv0Bits = randi([0 0],1,intv0);
    intv1Bits = randi([0 0],1,intv1);
    intv2Bits = randi([0 0],1,intv2);
    intv3Bits = randi([0 0],1,intv3);
    headBits  = randi([0 0],1,hLen);
    tailBits  = randi([0 0],1,tLen);
    payloadBits = randi([0 1],1,pLen);
    frame_hBits1 = [0 1 0 0 0 1 0 0 1 1 0 0 1 1 0 1] ;
    txBits1 = [frame_hBits1 intv0Bits frame_hBits1 intv1Bits  frame_hBits1 intv2Bits  frame_hBits1 intv3Bits payloadBits]
    txBits = [headBits txBits1 tailBits] ;
    txSymPr = mskmod(txBits,samplesPerSymbol);
    
    
    %txSym = [txSymPr,zeros(1,2048-192),txSymPr,zeros(1,2048-192),txSymPr,zeros(1,2048-192)]
    txSym=txSymPr;
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
        Frac = 10000; 
        timingOffset = round(ppm*Frac) ;
        rxSigTimingOff = resample(rxchan,Frac+timingOffset,Frac) ;
        %     varDelay = dsp.VariableFractionalDelay;
        %     rxSigTimingOff = varDelay(rxchan,0.2);
        %pwelch(rxchan,[],[],[],192e6,'twosided');  
        if flag_plot
            pwelch(rxSigTimingOff,[],[],[],192e6,'twosided');
            stem(real(rxchan(1:192)));
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
        rxSigCFO = rxSigTimingOff.*exp(-1*j*2*pi*freqOffset(ii)*[0:Len-1]*Ts) ;
    else
        rxSigCFO = rxSigTimingOff ;
    end    
    %
%     scatterplot(txSym(1:samplesPerSymbol:end));
%     scatterplot(rxSigCFO(1:samplesPerSymbol:end)) ;
    % Pass the signal through an AWGN channel
    
    if flag_shift
        rxSig = [rxSigCFO(shift+1:end),zeros(1,shift)];
    else
        rxSig = rxSigCFO ;
    end
    
    %%%get the corr for ppm estimation
    for tt=1:1
        rxCorr = rxSig((tt-1)*Ns*Nf+1:(tt)*Ns*Nf+0) ;
                
        tx1 = txSymPr(1:192) ;
%         plot(real(tx1),'g');
%         hold on; plot(imag(tx1),'b');
%         plot(angle(tx1),'r');
%         hold on;
%         plot(angle(rxCorr),'b');
    
        %freq_est = freq_est(tx1,Ns) ;
        f_est = freq_est(rxCorr',Ns) ;
        res_freq = (f_est*16e6-freqOffset(ii)) ;
        mse = (res_freq/16e6)^2 ; 
        fprintf("freq estimation is %f, ideal %f , error %f , res_offset %f angle %f mse %f \n" ,...
              f_est*16e6 , freqOffset(ii) , f_est*16e6/freqOffset(ii) , res_freq , res_freq *320/16e6*360 ,mse  ) ;
          
    end
    %
    % Save the transmitted signal for plotting
    plot_rx = rxSig;
    
%    
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	rx_demod=mskdemod(rxSig,12);
%     [err_num(ii) tt] = biterr(rx_demod(1:end-1),txBits(1:end-1)) ;
%     sum_err_num(ii) = sum_err_num(ii)+err_num(ii);
%     total_bits(ii) = total_bits(ii) + length(txBits) ;
%     
%     fprintf("SNR is %f and FRAME NUMBER is %d\n" , snr , jj ) ;
%     if( sum_err_num(ii) > 100 ) || (jj>1000) % at least 10-5
%         break ;
%     end  
end    
      
%     err_rate(ii) = sum_err_num(ii)/total_bits(ii) ;
%     ideal_rate(ii)=1/2*erfc(sqrt(10.^(snr/10)/2));
end

% if flag_plot

%     semilogy(SNR,err_rate);
%     hold on ;
%     semilogy(SNR,ideal_rate,'r');
% end
keyboard ;

