function   txsig =  signal_compensation(rxSig,freq,mphase,len);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	txsig=rxSig.*exp(-1*j*(2*pi*freq*[0:len-1]+mphase));
end
