function   location =  signal_detect(rxSigin,L);

max_papr = 0 ;
location = 0 ;


% Fpass = 0.1;             % Passband Frequency
% Fstop = 0.2;             % Stopband Frequency
% Dpass = 0.057501127785;  % Passband Ripple
% Dstop = 0.01;            % Stopband Attenuation
% dens  = 20;              % Density Factor
% 
% % Calculate the order from the parameters using FIRPMORD.
% [N, Fo, Ao, W] = firpmord([Fpass, Fstop], [1 0], [Dpass, Dstop]);
% 
% % Calculate the coefficients using the FIRPM function.
% b  = firpm(N, Fo, Ao, W, {dens});
% Hd = dfilt.dffir(b);
% 
% rxSig = filter(Hd.Numerator,1,rxSigin)
rxSig = rxSigin ;

papr_all = zeros(length(rxSig)-L+1,1) ;




for ii =1:length(rxSig)-L+1
    location = ii;
    curSig = rxSig(ii:ii+L-1);
    pxx=pwelch(curSig,[],[],[],192e6,'twosided');
    peak = mean([pxx(1:8); pxx(end-7:end)]) ;
    aver = mean(pxx(9:end-8)) ;
    papr = peak/aver;
    papr_all(ii) = papr ;
end
Numerator = ones(192,1)/192 ;
dd = filter(Numerator,1,papr_all) ;

[max_papr location] = max(dd) ;
location = location-96;
fprintf('location is: %d\n',location);
fprintf('papr is: %d\n',10*log(max_papr));

