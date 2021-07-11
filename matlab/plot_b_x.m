%Master Thesis R06323052 朝霍 

%plot 

%home production
t_sleep            = [516 533]; %何痸(と何) from Data_Time(new)
t_sleep_eat        = [588 596]; %何痸(と何)+ノ繺(甦) from Data_Time(new)
t_sleep_shower_eat = [644 662]; %何痸(と何)+縟瑍∟疍帝杆┪てЗ+ノ繺(甦)(ゲ璶丁) from Data_Time(new)
t_ness = t_sleep_shower_eat ; 

t_free_work = 24-t_ness/60 ;  % パ丁+丁
k0 = zeros(1,2) ; 
k1 = zeros(1,2) ; 
k2 = zeros(1,2) ; 
x0 = [187 164]; % from Data_Woman_new
x1 = [316 286]; % from Data_Woman_new
x2 = [375 332]; % from Data_Woman_new
x0 = x0/60./t_free_work;
x1 = x1/60./t_free_work;
x2 = x2/60./t_free_work;
for i = 1:2
  k0(i) = x0(i);
  k1(i) = x1(i)-k0(i);
  f=@(kk)(k0(i)+k1(i)*2^kk-x2(i));
  [kk, fval_k2]=fsolve(f,0.4845);
  k2(i) = kk;  
end

b = 0:0.01:3;
x_1 = k0(1)+k1(1).*b.^k2(1); 
x_2 = k0(2)+k1(2).*b.^k2(2); 
plot (b,x_1,'r-')
hold on;
plot (b,x_2,'b--')
hold on;
xlabel('b'), ylabel('x')
title('x=k_0+k_1*b^{k_2}')
hold on;
legend('1958~1962','1968~1972')