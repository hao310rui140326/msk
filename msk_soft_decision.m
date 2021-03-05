function  [freq,mphase] =  msk_soft_decision(rxSig,w,l,t,freq_all,phase_all);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for m=1:w
		txsig=signal_compensation(rxSig,freq_all(m),phase_all(m),l);
		txsig=msk_soft_decode(txsig,l,t);
		txsum(m)=(txsig*txsig')/(l-t+1);
	end
	txmax = max(txsum);
	idx = find(txsum==txmax);
	freq=freq_all(idx(1));	
	mphase=phase_all(idx(1));	
end
