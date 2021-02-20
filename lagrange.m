function  yi = lagrange(x,y,xi)

%  Lagrange 插值多项式，调用格式为

%   yi = LagInterp(x,y,xi)

%  其中

%  x 为插值节点，y为节点处函数值，

%  xi 为为估计函数自变量，yi 为xi处函数估计值

%

n = length(x); m = length(xi); p = zeros(n,m);

for k = 1:n

   t = ones(n,m);

   for j = 1:n

       if j~=k

           if abs(x(k) - x(j))<eps

              error('% 输入的插值节点必须互异！');

           end

           t(j,:) = (xi - x(j))/(x(k) - x(j));

       end

   p(k,:) = prod(t);

   end

   yi = y*p;

end