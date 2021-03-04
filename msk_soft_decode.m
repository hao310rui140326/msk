function   x =  msk_soft_decode(rxSig,len,t);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	rx_real = real(rxSig);
	rx_imag = imag(rxSig);
	for a=1:len-t+1
		b = 0:t-1;
		x_real(a)=sum(rx_real(a+b).*cos(b*pi/12));
		x_imag(a)=sum(rx_imag(a+b).*sin(b*pi/12));
	end
	x=complex(x_real,x_imag);
end
