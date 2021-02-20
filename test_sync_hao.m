clear all;
close all;
clc;

numFrames = 1;

samplesPerSymbol= 12;

symbolsPerFrame = 16 ;

fSample = 192e6;
fSymbol = 16e6 ;

fs = fSample ;
Ts = 1.0/fs ;

timingOffset = 0;
varDelay = dsp.VariableFractionalDelay;

freqOffset = 5e3;
phaseOffset = 0;
pfo = comm.PhaseFrequencyOffset(...
    'FrequencyOffset',freqOffset, ...
    'PhaseOffset',phaseOffset, ...
    'SampleRate',samplesPerSymbol/Ts);

%EbNo = SNR + 10*log10(samplesPerSymbol);
EbNo = 0 ; 

chAWGN = comm.AWGNChannel(...
    'NoiseMethod','Signal to noise ratio (Eb/No)', ...
    'EbNo',EbNo,...
    'SignalPower',1, ...
    'SamplesPerSymbol',samplesPerSymbol);

timeSync = comm.MSKTimingSynchronizer(...
    'SamplesPerSymbol',samplesPerSymbol, ...
    'ErrorUpdateGain',0.02);

phaseSync = comm.CarrierSynchronizer(...
    'Modulation','QPSK', ...
    'ModulationPhaseOffset','Custom', ...
    'CustomPhaseOffset',0, ...
    'SamplesPerSymbol',1);


    %------------------------------------------------------------------------
    % Generate and modulate data
    %------------------------------------------------------------------------
    %txBits = randi([0 1],symbolsPerFrame,1);
    headBits = randi([0 1],1,92);
    tailBits = randi([0 1],1,92);
    txBits1 = [0 1 0 0 0 1 0 0 1 1 0 0 1 1 0 1] ;
    txBits = [headBits txBits1 tailBits] ;
    txSym = mskmod(txBits,samplesPerSymbol);
    %------------------------------------------------------------------------
    % Transmit through channel
    %------------------------------------------------------------------------
    %
    % Add timing offset
    rxSigTimingOff = varDelay(txSym,timingOffset*samplesPerSymbol);
    %
    % Add carrier frequency and phase offset
    rxSigCFO = pfo(rxSigTimingOff);
    %
    scatterplot(txSym(1:samplesPerSymbol:end));
    scatterplot(rxSigCFO(1:samplesPerSymbol:end)) ;
    % Pass the signal through an AWGN channel
    rxSig = chAWGN(rxSigCFO);
    %
    % Save the transmitted signal for plotting
    plot_rx = rxSig;
    
    plot(abs(fft(rxSig.^4,256)));
    
    rx_demod=mskdemod(rxSig,12);
    figure;plot(rx_demod(93:93+15)-txBits1);
   

    for i=1:(length(rx_demod)-16)
        %%rx_mult(i)=(rx_demod(i:i+15)*txBits1');
        %%rx_mult(i)=xcorr(rx_demod(i:i+15),txBits1);
        rx_mult(i)=0;
        for j=1:16
            if (rx_demod(i+j-1)==txBits1(j)) 
                rx_mult(i)=rx_mult(i)+1;
            else
                rx_mult(i)=rx_mult(i)-1;
            end
        end
    end
    figure;plot(rx_mult);


    %
% %     %------------------------------------------------------------------------
% %     % Timing recovery
% %     %------------------------------------------------------------------------
% %   recoverTimingPhase = 1;
% %     if recoverTimingPhase
% %         % Recover symbol timing phase using fourth-order nonlinearity
% %         % method
% %         [rxSym,timEst] = timeSync(rxSig);
% %         % Calculate the timing delay estimate for each sample
% %         timEst = timEst(1)/samplesPerSymbol;
% %     else
% %         % Do not apply timing recovery and simply downsample the received
% %         % signal
% %         rxSym = downsample(rxSig,samplesPerSymbol);
% %         timEst = 0;
% %     end
% % 
% %     % Save the timing synchronized received signal for plotting
% %     plot_rxTimeSync = rxSym;
% % 
% %     %------------------------------------------------------------------------
% %     % Carrier frequency and phase recovery
% %     %------------------------------------------------------------------------
% %     recoverCarrier = 1;
% %     if recoverCarrier
% %         % The following script applies carrier frequency and phase recovery
% %         % using a second order phase-locked loop (PLL), and removes phase ambiguity
% %         [rxSym,phEst] = phaseSync(rxSym);
% %         removePhaseAmbiguityMSKSignalRecoveryEx;
% %         freqShiftEst = mean(diff(phEst)/(Ts*2*pi));
% %         phEst = mod(mean(phEst),360); % in degrees
% %     else
% %         freqShiftEst = 0;
% %         phEst = 0;
% %     end
% % 
% %     % Save the phase synchronized received signal for plotting
% %     plot_rxPhSync = rxSym;
% %     %------------------------------------------------------------------------
% %     % Demodulate the received symbols
% %     %------------------------------------------------------------------------
% %     rxBits = mskdemod(rxSym,1);
% %     %------------------------------------------------------------------------
% %     % Calculate the bit error rate
% %     %------------------------------------------------------------------------
% %     errorStats = BERCalc(txBits,rxBits);
% %     %------------------------------------------------------------------------
% %     % Plot results
% %     %------------------------------------------------------------------------
% %     plotResultsMSKSignalRecoveryEx;
keyboard ;

