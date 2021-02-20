function t_est = time_est(data,N) ;


rm1 = zeros(N,1) ;
rm2 = zeros(N,1) ;
rm3 = zeros(N,1) ;
rm4 = zeros(N,1) ;
rm5 = zeros(N,1);
rm6 = zeros(N,1);
rm7 = zeros(N,1);
rm8 = zeros(N,1);

rm1 = (rm_calc(data,N,1)) ;
%rm2 = (rm_calc(data,N,2)) ;
%rm3 = (rm_calc(data,N,3)) ;
%rm4 = (rm_calc(data,N,4)) ;
%rm5 = (rm_calc(data,N,5)) ;
%rm6 = (rm_calc(data,N,6)) ;
%rm7 = (rm_calc(data,N,7)) ;
%rm8 = (rm_calc(data,N,8)) ;







sum_res_all = 0 ;
for ii = 0:N-1
    rm_aver = rm1(ii+1) + rm2(ii+1) + rm3(ii+1) + rm4(ii+1) ...
        + rm5(ii+1) + rm6(ii+1) + rm7(ii+1) + rm8(ii+1) ;
    sum_res_all = sum_res_all + abs(rm_aver) * exp(-1*j*2*pi*(ii)/N) ;
end    

theta = angle(sum_res_all) ; 

t_est = -1  * theta/(2*pi);   
