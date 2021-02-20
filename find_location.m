function   location =  find_location(rxSig,txBits1,err);
%%find location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%disp('location is :');
	rx_demod(1,:)=mskdemod(rxSig,12);
    for m=1:11
    	rx_demod(m+1,:)=mskdemod([rxSig(m:1:end),zeros(1,m-1)],12);
    end
    
    location=[];
    
    for i=1:12
    	for m=1:(length(rx_demod(1,:))-15)
        	%%rx_mult(i)=(rx_demod(i:i+15)*txBits1');
        	%%rx_mult(i)=xcorr(rx_demod(i:i+15),txBits1);
        	rx_mult(i,m)=0;
        	for j=1:16
        	    if (rx_demod(i,m+j-1)==txBits1(j)) 
        	        rx_mult(i,m)=rx_mult(i,m)+1;
        	    else
        	        rx_mult(i,m)=rx_mult(i,m)-1;
        	    end
            end
        end
    	%%subplot(1,12,i);plot(rx_mult(i,:));
        nlocation = find(rx_mult(i,:)>=err);
        %%disp(nlocation);
        for n =1:length(nlocation)
            idx=find(location==nlocation(n));
            if idx>0 
                location = location;
            else
                location = [location,(nlocation(n)-1)*12+i];
            end
        end
    end

    location=sort(location);
   disp('location is: ');
   disp(location); 
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




end
