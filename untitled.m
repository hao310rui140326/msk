
ppm_a = [0:0.005:0.25] ;
for ii = 1:length(ppm_a)
    ppm = ppm_a(ii) ;
adata = [0:pi/8:pi-pi/8]+ppm*pi/2*[0:7] ;

adata0 = adata(1:4) ;
adata1 = adata(5:8) ;

data = cos(adata) + j*sin(adata) ; 
data0 = cos(adata0) + j*sin(adata0) ;
data1 = cos(adata1) + j*sin(adata1) ;

R = (conj(data0).*(data1)).^2 ;



res =  0 ;
for ii = 1:4
    res = res + R(ii) ;% * exp(-j*2*pi*(ii-1)/4) ;
end

ppm_est = angle(res*-1/4)/4/pi;

fprintf("ideal ppm is %f , real ppm is %f\n" , ppm , ppm_est) ;

end
keyboard;



