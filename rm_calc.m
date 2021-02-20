function rm_calc_res = rm_calc(data,N,m) ;

len = length(data) ;
L0 = len/N ;
%% m must less than L0

sum_res = zeros(N,1) ;
for ii = 0: N-1
    for jj = m+1: L0
        xk  = data(ii+1+N*(jj-1)) ;
        xkm = conj( data(ii+1+N*(jj-1-m)) ) ;
        sum_res(ii+1) = sum_res(ii+1) + ((xk*xkm)^2)/(L0-m) ;
        %fprintf("ii %d and jj %d add item is %f\n",ii,jj,((xk*xk1)^2)) ;
    end    
    %fprintf("ii %d est abs %f is %f + j * %f \n",ii,abs(sum_res(ii+1)) ,real(sum_res(ii+1)) , imag(sum_res(ii+1))) ;
end    



rm_calc_res = sum_res ; 
