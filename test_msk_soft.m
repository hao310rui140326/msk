rxSig = 1:192;
w=10;
l=192;
t=12;
freq_all=1:10;
phase_all=1:10;


[freq,phase]=msk_soft_decision(rxSig,w,l,t,freq_all,phase_all);


