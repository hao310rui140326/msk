function f_est = freq_est(data,N) ;

len = length(data) ;
L0 = len/N ;

sum_res_all = 0;

sum_res = zeros(N,1) ;
for ii = 0: N-1
    for jj = 1: L0-1 
        xk  = data(ii+1+N*(jj)) ;
        xk1 = conj( data(ii+1+N*(jj-1)) ) ;
        sum_res(ii+1) = sum_res(ii+1) + ((xk*xk1)^2)/(L0-1) ;
    end    
end 

sum_res2 = zeros(N,1) ;
for ii = 0: N-1
    for jj = 0: L0-1 
        xk  = data(ii+1+N*(jj)) ;
        xk1 = conj( data(ii+1+N*(jj)) ) ;
        sum_res2(ii+1) = sum_res2(ii+1) + ((xk*xk1)^2)/(L0) ;
    end    
end  


%%for ii=1:N
%%    theta(ii) = angle( -sum_res(ii) * conj(sum_res2(ii)) ) ;  
%%    f_est(ii) =    theta(ii)/(4*pi);   
%%end

    theta = angle( -sum_res(1) * conj(sum_res2(1)) ) ;  
    f_est =    theta/(4*pi);   
