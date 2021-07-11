%Master Thesis R06323052 陳思齊 
%calibration (wage deifinition is only wage(excluding owner))
clear
clc

%price parameter


%price parameter

%1 definition(受雇人員報酬)
am_1      = [1.0000 1.1700 1.2835];  %from calibration 價格相關變數(受雇人員報酬)
wf_L_1    = [0.5235 0.6073 0.6806];  %from calibration 價格相關變數(受雇人員報酬)
wf_H_1    = [1.0031 1.0999 1.1490];  %from calibration 價格相關變數(受雇人員報酬)
%original(including 1822) 
pb_1_o    = [0.0912 0.1327 0.1639]; 
pe_1_o    = [0.0628 0.0632 0.0559];
c_t_1_ao  = [0.6467 0.7444 0.8335]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(total)
c_L_1_ao  = [0.6045 0.6795 0.7494]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(高中職以下)
c_H_1_ao  = [0.8908 0.9781 1.0306]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(大專以上)
c_t_1_to  = [0.7333 0.8232 0.8966]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(total)
c_L_1_to  = [0.6524 0.7126 0.7618]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(高中職以下)
c_H_1_to  = [1.0941 1.1444 1.1419]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(大專以上)
c_t_1_aom = [0.7711 0.8729 0.9589]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(total)
c_L_1_aom = [0.7288 0.8083 0.8771]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(高中職以下)
c_H_1_aom = [1.0536 1.1340 1.1504]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(大專以上)
c_t_1_tom = [0.8593 0.9429 1.0070]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(total)
c_L_1_tom = [0.7793 0.8382 0.8858]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(高中職以下)
c_H_1_tom = [1.2155 1.2470 1.2275]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(大專以上)
%new(excluding 1822)
pb_1_n   = [0.0829 0.0981 0.1523];
pe_1_an  = [0.0911 0.1021 0.1094]; % a: average(分層加權)
c_t_1_an = [0.6315 0.7368 0.8161]; % a: average(分層加權)(total)
c_L_1_an = [0.5920 0.6758 0.7325]; % a: average(分層加權)(高中職以下)
c_H_1_an = [0.8590 0.9571 1.0078]; % a: average(分層加權)(大專以上)
pe_1_tn  = [0.0984 0.1089 0.1152]; % t: total(全部平均)
c_t_1_tn = [0.7117 0.8125 0.8804]; % t: total(全部平均)(total)
c_L_1_tn = [0.6361 0.7075 0.7445]; % t: total(全部平均)(高中職以下)
c_H_1_tn = [1.0484 1.1172 1.1277]; % t: total(全部平均)(大專以上)

%2 definition(受雇人員報酬+產業主所得)
am_2      = [1.0000 1.1681 1.2772];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
wf_L_2    = [0.4849 0.5751 0.6391];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
wf_H_2    = [0.9242 1.0233 1.0671];  %from calibration 價格相關變數(受雇人員報酬+產業主所得)
%original(including 1822) 
pb_2_o    = [0.0835 0.1187 0.1494]; 
pe_2_o    = [0.0578 0.0589 0.0539];
c_t_2_ao  = [0.6045 0.6939 0.7686]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(total)
c_L_2_ao  = [0.5723 0.6432 0.6997]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(高中職以下)
c_H_2_ao  = [0.8209 0.9038 0.9439]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(大專以上)
c_t_2_to  = [0.6748 0.7554 0.8182]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(total)
c_L_2_to  = [0.6141 0.6711 0.7104]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(高中職以下)
c_H_2_to  = [1.0070 1.0492 1.0469]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(大專以上)
c_t_2_aom = [0.7221 0.8129 0.8836]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(total)
c_L_2_aom = [0.6901 0.7624 0.8169]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(高中職以下)
c_H_2_aom = [0.9739 1.0444 1.0524]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(大專以上)
c_t_2_tom = [0.7929 0.8671 0.9210]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(total)
c_L_2_tom = [0.7334 0.7872 0.8243]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(高中職以下)
c_H_2_tom = [1.1185 1.1452 1.1263]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(大專以上)

%newexcluding 1822)
pb_2_n   = [0.0737 0.1035 0.1364];
pe_2_an  = [0.0852 0.0948 0.1016]; % a: average(分層加權)
c_t_2_an = [0.5934 0.6846 0.7581]; % a: average(分層加權)(total)
c_L_2_an = [0.5634 0.6368 0.6896]; % a: average(分層加權)(高中職以下)
c_H_2_an = [0.7936 0.8831 0.9283]; % a: average(分層加權)(大專以上)
pe_2_tn  = [0.0911 0.1002 0.1062]; % t: total(全部平均)
c_t_2_tn = [0.6583 0.7436 0.8085]; % t: total(全部平均)(total)
c_L_2_tn = [0.6017 0.6635 0.6999]; % t: total(全部平均)(高中職以下)
c_H_2_tn = [0.9682 1.0223 1.0391]; % t: total(全部平均)(大專以上)

%3 definition(經常性收入)
am_3      = [1.0000 1.1773 1.2917];  %from calibration 價格相關變數(經常性收入)
wf_L_3    = [0.4612 0.5550 0.6219];  %from calibration 價格相關變數(經常性收入)
wf_H_3    = [0.8671 0.9702 1.0157];  %from calibration 價格相關變數(經常性收入)
%original(including 1822) 
pb_3_o    = [0.0724 0.1028 0.1295]; 
pe_3_o    = [0.0501 0.0510 0.0467];
c_t_3_ao  = [0.5235 0.6012 0.6659]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(total)
c_L_3_ao  = [0.4957 0.5573 0.6062]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(高中職以下)
c_H_3_ao  = [0.7111 0.7831 0.8179]; % a: average(分層加權) ((生活開銷-子女部分)/成年人*2)(大專以上)
c_t_3_to  = [0.5844 0.6546 0.7089]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(total)
c_L_3_to  = [0.5318 0.5815 0.6155]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(高中職以下)
c_H_3_to  = [0.8723 0.9091 0.9071]; % t: total(全部平均)   ((生活開銷-子女部分)/成年人*2)(大專以上)
c_t_3_aom = [0.6254 0.7043 0.7657]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(total)
c_L_3_aom = [0.5976 0.6605 0.7078]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(高中職以下)
c_H_3_aom = [0.8437 0.9049 0.9119]; % a: average(分層加權) (僅考慮子女(<23歲)+夫妻家庭)(大專以上)
c_t_3_tom = [0.6867 0.7513 0.7981]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(total)
c_L_3_tom = [0.6351 0.6821 0.7142]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(高中職以下)
c_H_3_tom = [0.9689 0.9923 0.9760]; % t: total(全部平均)   (僅考慮子女(<23歲)+夫妻家庭)(大專以上)
%new(excluding 1822)
pb_3_n   = [0.0639 0.0898 0.1184];
pe_3_an  = [0.0738 0.0822 0.0880]; % a: average(分層加權)
c_t_3_an = [0.5139 0.5929 0.6564]; % a: average(分層加權)(total)
c_L_3_an = [0.4879 0.5516 0.5970]; % a: average(分層加權)(高中職以下)
c_H_3_an = [0.6875 0.7650 0.8041]; % a: average(分層加權)(大專以上)
pe_3_tn  = [0.0790 0.0868 0.0920]; % t: total(全部平均)
c_t_3_tn = [0.5702 0.6441 0.7001]; % t: total(全部平均)(total)
c_L_3_tn = [0.5212 0.5748 0.6060]; % t: total(全部平均)(高中職以下)
c_H_3_tn = [0.8387 0.8856 0.9000]; % t: total(全部平均)(大專以上)


%control part
am   =     am_1 ; 
wf_L =   wf_L_1 ; 
wf_H =   wf_H_1 ;
pb   =   pb_1_o ;     
pe   =   pe_1_o ;    
c_L  = c_L_1_tn ;    
c_H  = c_H_1_tn ;


%home production
t_sleep            = [516 524 533]; %睡眠(含午睡) from Data_Time(new)
t_sleep_shower     = [565 572 574]; %睡眠(含午睡)+盥洗、沐浴、著裝 from Data_Time(new)
t_sleep_eat        = [588 596 596]; %睡眠(含午睡)+用餐(含吃宵夜) from Data_Time(new)
t_sleep_shower_eat = [644 653 662]; %睡眠(含午睡)+盥洗、沐浴、著裝或化妝+用餐(含吃宵夜)(必要時間) from Data_Time(new)
t_ness = t_sleep ; 

t_free_work = 24-t_ness/60 ;  % 自由時間+工作時間
k0 = zeros(1,3) ; 
k1 = zeros(1,3) ; 
k2 = zeros(1,3) ; 
x0 = [187 178 164]; % from Data_Woman_new
x1 = [316 308 286]; % from Data_Woman_new
x2 = [375 354 332]; % from Data_Woman_new
x0 = x0/60./t_free_work;
x1 = x1/60./t_free_work;
x2 = x2/60./t_free_work;
for i = 1:3
  k0(i) = x0(i);
  k1(i) = x1(i)-k0(i);
  f=@(kk)(k0(i)+k1(i)*2^kk-x2(i));
  [kk, fval_k2]=fsolve(f,0.4845);
  k2(i) = kk;  
end


%preference parameter

%fertility related
b_L      = [2.253 2.053 1.748];  %from 各世代平均生育總數(66)
b_H      = [1.465 1.402 1.402];  %from 各世代平均生育總數(66)
b_total  = [2.213 2.003 1.781];  %from 各世代平均生育總數(66)
R        = [0.435 0.432 0.431];  %from PSFD (L->L percentage)

%production side parameter
tf_L = [143 170 190]; %from Data_Woman_work (min./day)
tf_H = [236 239 255]; %from Data_Woman_work (min./day)
tf_L = tf_L/60./t_free_work ; %標準化(換算成小時後，在換算成自由時間占比)
tf_H = tf_H/60./t_free_work ; %標準化(換算成小時後，在換算成自由時間占比)
n_L  = [5823989 5495228 4702371]; % from 各世代平均生育總數
n_H  = [ 809586 1143609 1691866]; % from 各世代平均生育總數
N_L  = n_L./(n_L+n_H) ; %標準化N_L
N_H  = n_H./(n_L+n_H) ; %標準化N_H
Z_L  = N_L.*tf_L ; %L_type labor force
Z_H  = N_H.*tf_H ; %H_type labor force

% Assume
rho = 0.401 ; % from Krusell(2000)



%calibration for gamma phi eta  & part of beta A pi
bd_L_m      = zeros(1,3);
be_L_m      = zeros(1,3);
b_L_m       = zeros(1,3);
R_m         = zeros(1,3);
b_total_m   = zeros(1,3);
tfd_L_m     = zeros(1,3);
tfe_L_m     = zeros(1,3);
tf_L_m      = zeros(1,3);
tf_H_m      = zeros(1,3);
Z_L_m       = zeros(1,3);
Z_H_m       = zeros(1,3);
beta        = zeros(1,3);
A           = zeros(1,3);
pi          = zeros(1,3);

for i = 1
x0_t = [b_L(i) b_L(i) b_H(i) R(i) 1.1418 1.3558 0.7686 1.1885 0.5454]; % bd_L be_L b_H R gamma phi eta A pi
N = 100000 ;
options = optimset('MaxFunEvals', N, 'MaxIter', N, 'Algorithm', 'levenberg-marquardt');
[x,fval_cali]=fsolve(@(x)myfun_2at(x,pb(i),pe(i),am(i),c_L(i),c_H(i),k0(i),k1(i),k2(i),wf_L(i),wf_H(i),tf_L(i),tf_H(i),b_L(i),b_H(i),N_L(i),N_H(i),rho),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    R_m(i)       = x(3) ;
    b_L_m(i)     = R_m(i)*bd_L_m(i)+(1-R_m(i))*be_L_m(i) ;
    b_total_m(i) = N_L(i)*b_L_m(i)+N_H(i)*b_H(i) ;
    gamma        = x(4) ;
    alpha_L      = x(5) ;
    alpha_H      = x(6) ;
    A(i)         = x(7) ; 
    pi(i)        = x(8) ; 
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb(i)+pe(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb(i)+pe(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
    beta(i)      = (-((1+gamma)*log((pb(i)*bd_L_m(i)+wf_L(i)*k1(i)*k2(i)*(bd_L_m(i))^(k2(i)))/((pb(i)+pe(i))*be_L_m(i)+wf_L(i)*k1(i)*k2(i)*(be_L_m(i))^(k2(i))))+alpha_L*log(bd_L_m(i)/be_L_m(i)))/(log(wf_L(i)/wf_H(i)))) ;
calibration_P1    = [alpha_L alpha_H gamma];
disp('      alp_L      alp_H   gamma')
disp(calibration_P1)
calibration_P2    = [beta(i) A(i) pi(i) R_m(i)];
disp('      beta          A             pi          R_m')
disp(calibration_P2)
end



%caliabration for part of beta A pi & part of pb pe c_
bd_L_m      = zeros(1,3);
be_L_m      = zeros(1,3);
b_L_m       = zeros(1,3);
R_m         = zeros(1,3);
b_total_m   = zeros(1,3);
tfd_L_m     = zeros(1,3);
tfe_L_m     = zeros(1,3);
tf_L_m      = zeros(1,3);
tf_H_m      = zeros(1,3);
Z_L_m       = zeros(1,3);
Z_H_m       = zeros(1,3);

beta   = zeros(1,3);
A      = zeros(1,3);
pi     = zeros(1,3);
pb_m   = zeros(1,3);
pe_m   = zeros(1,3);
c_L_m  = zeros(1,3);
c_H_m  = zeros(1,3);

%%
for i = 1:3
x0_t = [b_L(i) b_L(i) 0.4 pb(i) pe(i) c_L(i) c_H(i) 1.1885 0.5454 ];  
N = 100000 ;
options = optimset('MaxFunEvals', N, 'MaxIter', N, 'Algorithm', 'levenberg-marquardt');
[x,fval_cali]=fsolve(@(x)myfun_2a2ct(x,alpha_L,alpha_H,am(i),gamma,k0(i),k1(i),k2(i),wf_L(i),wf_H(i),b_L(i),b_H(i),N_L(i),N_H(i),rho,tf_L(i),tf_H(i),beta(i)),x0_t,options);
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    R_m(i)       = x(3) ;
    b_L_m(i)     = R_m(i)*bd_L_m(i)+(1-R_m(i))*be_L_m(i) ;
    b_total_m(i) = N_L(i)*(b_L_m(i))+N_H(i)*b_H(i) ;
    pb_m(i)      = x(4) ;
    pe_m(i)      = x(5) ;
    c_L_m(i)     = x(6) ;
    c_H_m(i)     = x(7) ;
    A(i)         = x(8) ; 
    pi(i)        = x(9) ; 
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb_m(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
    calibration_P2    = [beta(i) A(i) pi(i) R_m(i)];
    disp('      beta          A             pi          R_m')
    disp(calibration_P2)
end
disp('     pb_m     pb_m      pb_m')
disp(pb_m)
disp('     pe_m     pe_m      pe_m')
disp(pe_m)
disp('    c_L_m      c_L_m       c_L_m')
disp(c_L_m)

%%
%caliabration_check
bd_L_m      = zeros(1,3);
be_L_m      = zeros(1,3);
b_L_m       = zeros(1,3);
b_H_m       = zeros(1,3);
b_total_m   = zeros(1,3);
wf_L_m      = zeros(1,3);
wf_H_m      = zeros(1,3);
R_m         = zeros(1,3);
tfd_L_m     = zeros(1,3);
tfe_L_m     = zeros(1,3);
tf_L_m      = zeros(1,3);
tf_H_m      = zeros(1,3);
lfd_L_m     = zeros(1,3);
lfe_L_m     = zeros(1,3);
lf_L_m      = zeros(1,3);
lf_H_m      = zeros(1,3);
cd_L_m      = zeros(1,3);
ce_L_m      = zeros(1,3);
c_L_m       = zeros(1,3);
c_H_m       = zeros(1,3);
Z_L_m       = zeros(1,3);
Z_H_m       = zeros(1,3);

%pb_m = 1.01*pb_m  ;
%k1 = k1*0.985 ;
   Q = 0.0 ;
for i = 1:3
    a = i ; 
    b = i ;
    c = i ;
    d = i ;
    e = i ;
xi = [b_L(i), b_L(i), b_H(i), R(i), wf_L(i), wf_H(i)] ;
N = 1000000 ;
options = optimset('MaxFunEvals',N,'MaxIter',N, 'Algorithm', 'levenberg-marquardt');
% j == >  preference parameter 選取變數(for phi eta gamma)
[x,fval_m]=fsolve(@(x)myfun_2a_c(x,pb_m(a),pe_m(b),am(c),c_m(d),k0(e),k1(e),k2(e),N_L(i),N_H(i),rho,alpha_L,alpha_H,gamma,beta(i),A(i),pi(i)),xi,options);
    pb_m(i)      = pb_m(a)*(1-Q) ;
    pe_m(i)      = pe_m(b) ;
    am(i)        = am(c)-(N_L(i)*(x(4)*x(1)+(1-x(4))*x(2)))*pb_m(i)*Q ;
    c_m(i)       = c_m(d) ;
    k0(i)        = k0(e) ;
    k1(i)        = k1(e) ;
    k2(i)        = k2(e) ;
    bd_L_m(i)    = x(1) ;
    be_L_m(i)    = x(2) ;
    b_L_m(i)     = (x(4)*x(1)+(1-x(4))*x(2)) ;
    b_H_m(i)     = x(3) ;
    R_m(i)       = x(4) ;
    b_total_m(i) = (N_L(i)*(x(4)*x(1)+(1-x(4))*x(2))+N_H(i)*x(3)) ;
    wf_L_m(i)    = x(5) ;
    wf_H_m(i)    = x(6) ;
    tfd_L_m(i)   = (1-k0(i)-((gamma*pb_m(i)*bd_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tfe_L_m(i)   = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*be_L_m(i)+(gamma*k2(i)+alpha_L)*wf_L(i)*k1(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L(i)))) ;
    tf_L_m(i)    = R_m(i)*tfd_L_m(i)+(1-R_m(i))*tfe_L_m(i) ;
    tf_H_m(i)    = (1-k0(i)-((gamma*(pb_m(i)+pe_m(i))*b_H(i)+(gamma*k2(i)+alpha_H)*wf_H(i)*k1(i)*b_H(i)^k2(i))/(alpha_H*wf_H(i)))) ;
    lfd_L_m(i)   = 1-tfd_L_m(i)-k0(i)-k1(i)*bd_L_m(i)^k2(i) ;
    lfe_L_m(i)   = 1-tfe_L_m(i)-k0(i)-k1(i)*be_L_m(i)^k2(i) ;
    lf_L_m(i)    = R_m(i)*lfd_L_m(i)+(1-R_m(i))*lfe_L_m(i) ;
    lf_H_m(i)    = 1-tf_H_m(i)-k0(i)-k1(i)*b_H_m(i)^k2(i) ;
    cd_L_m(i)    = (pb_m(i)*bd_L_m(i)+wf_L_m(i)*k1(i)*k2(i)*bd_L_m(i)^k2(i))/(alpha_L*wf_L_m(i)) ;
    ce_L_m(i)    = ((pb_m(i)+pe_m(i))*be_L_m(i)+wf_L_m(i)*k1(i)*k2(i)*be_L_m(i)^k2(i))/(alpha_L*wf_L_m(i)) ;
    c_L_m(i)     = R_m(i)*cd_L_m(i)+(1-R_m(i))*ce_L_m(i) ;
    c_H_m(i)    = ((pb_m(i)+pe_m(i))*b_H_m(i)+wf_H_m(i)*k1(i)*k2(i)*b_H_m(i)^k2(i))/(alpha_H*wf_H_m(i)) ;
    Z_L_m(i)     = N_L(i)*tf_L_m(i) ;
    Z_H_m(i)     = N_H(i)*tf_H_m(i) ;
end
disp('      R_m        R_m       R_m')
disp(R_m)
disp('       b_L          b_L         b_L')
disp(b_L)
%disp('   bd_L_m bd_L_m  bd_L_m')
%disp(bd_L_m)
%disp('   be_L_m  be_L_m  be_L_m')
%disp(be_L_m)
disp('    b_L_m    b_L_m    b_L_m')
disp(b_L_m)
disp('       b_H         b_H        b_H')
disp(b_H)
disp('   b_H_m    b_H_m   b_H_m')
disp(b_H_m)
%disp('       pb            pb           pb')
%disp(pb)
disp('     pb_m      pb_m     pb_m')
disp(pb_m)
%disp('        pe           pe           pe')
%disp(pe)
disp('     pe_m      pe_m     pe_m')
disp(pe_m)
%disp('        c_            c_            c_')
%disp(c_)
%%
disp('      c_m         c_m        c_m')
disp(c_m)
disp('    b_t_m     b_t_m     b_t_m')
disp(b_total_m)
disp('      lfd_L        lfd_L      lfd_L')
disp(lfd_L_m)
disp('     tfd_L        tfd_L       tfd_L')
disp(tfd_L_m)
disp('      cd_L        cd_L      cd_L')
disp(cd_L_m)
disp('      lfe_L        lfe_L      lfe_L')
disp(lfe_L_m)
disp('     tfe_L        tfe_L       tfe_L')
disp(tfe_L_m)
%disp('       tf_L         tf_L          tf_L')
disp('      ce_L        ce_L        ce_L')
disp(ce_L_m)
%disp(tf_L)
disp('    lf_L_m    lf_L_m     lf_L_m')
disp(lf_L_m)
disp('    tf_L_m    tf_L_m    tf_L_m')
disp(tf_L_m)
disp('    c_L_m     c_L_m    c_L_m')
disp(c_L_m)
%disp('       tf_H         tf_H          tf_H')
%disp(tf_H)
disp('   lf_H_m    lf_H_m   lf_H_m')
disp(lf_H_m)
disp('   tf_H_m    tf_H_m   tf_H_m')
disp(tf_H_m)
disp('   c_H_m    c_H_m    c_H_m')
disp(c_H_m)
disp('   wf_L_m  wf_L_m   wf_L_m')
disp(wf_L_m)
disp('   wf_H_m wf_H_m wf_H_m')
disp(wf_H_m)
disp('    am_m     am_m     am_m')
disp(am)