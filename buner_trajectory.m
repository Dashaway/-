%181016
%Dashaway
%燃烧参数计算
%圆柱装药，等面燃烧

%190117
%增加注释
%换为圆柱药柱，药柱数量增加

%201101
%换为单管多室发动机

%201211
%增加主动段外弹道

%201229
%修复了一些bug，增加注释

%210108
%多次点火

%210113
%扫描循环

clear;
close all;
%初始值设定
FileName = 'C_xon.xlsx';
[data] = xlsread(FileName,'Sheet1', 'B1:D165' );
NUM = xlsread(FileName,'Sheet1', 'A2:A2' );
%扫描循环设定
%角度循环设定
theta_sta = 30;       %初始射角(deg)
theta_ste = 5;       %步进射角(deg)
theta_end = 75;       %结束射角(deg)
%点火循环设定
t_sta = 0;      %续航药柱点火初始时间(s)
t_ste = 0.5;      %续航药柱点火步进时间(s)
t_end = 30;      %续航药柱点火结束时间(s)

%常量

%计算辅助常量
long = 1000000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 1e-4;      %步进值
%燃烧室参数
Dr_t = 0.03;       %推进药柱燃烧室内径(m)
Dr_x = 0.03;       %续航药柱燃烧室内径(m)
At = 9.62e-4;     %喷管喉部面积(m^2)
D_p = 0.12;      %弹径(m)
m_k = 7;         %起飞弹重(kg)
Vg_z = 0;       %未装药空间体积(m^3)
%装药参数
D_t = 0.03;         %推进药柱初始外径(m)
d_t = 0.01;       %推进药柱初始内径(m)
D_x = 0.03;         %续航药柱初始外径(m)
d_x = 0.01;       %续航药柱初始内径(m)
n_t = 2;        %推进药柱数量
n_x = 2;        %续航药柱数量
Lp = 0.25;      %装药长度(m)
%已知常量
p0 = 101325;    %初始压强(Pa)
gamma = 1.2;        %比热比
rho_p = 1730;       %密度(kg/m^3)
n_p = 0.302;      %压强指数
c = 1600;       %特征速度(m/s)
rb_0 = 5e-3;      %?初始燃速(m/s)
alpha_r = 1.7e-4;        %?燃速系数
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?
g = 9.78030;      %重力加速度(m/s^2)
R_1 = 287;        %空气气体常数(J/kg/K)
rho_y0 = 1.293;       %空气密度(kg/m^3)
v_e = 720;       %燃气出口速度(m/s)
%外弹道参数
v_0 = 1;       %出口速度(m/s)
x_0 = 0;           %水平初始位置(m)
y_0 = 1.1;           %竖直初始位置(m)
tau_0 = 289.1;        %初始虚温(K) 
y_end = 400;           %竖直结束位置(m)

%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
s_t0 = n_t*pi*d_t;      %推进药柱燃烧面初始周长(m)
s_x0 = n_x*pi*d_x;      %续航药柱燃烧面初始周长(m)
Ap_t0 =  n_t*pi*d_t^2 / 4;        %推进药柱初始通气面积(m^2)
Vp_t0 = (n_t*pi*(D_t^2) / 4 - Ap_t0)*Lp;       %推进药柱药柱总体积(m^3)
Vp_x0 = (n_x*pi*(D_x^2) / 4 - n_x*pi*d_x^2 / 4)*Lp;       %续航药柱药柱总体积(m^3)
mp = (Vp_t0 + Vp_x0)*rho_p;     %药柱总质量(kg)
c_d0 = 1000*D_p*D_p/m_k;      %初始弹道系数
c_s0 = (1.4*R_1*tau_0)^0.5;          %初始声速(m/s)

%约束条件

%推进段结束条件
ep_t = (D_t - d_t)/2;     %推进段肉厚(m)

%续航段结束条件
ep_x = (D_x - d_x)/2;     %续航段肉厚(m)


%变量
%计算用变量（已赋初值）
p = p0*ones(1,long);        %实际压强(Pa)
r = rb_0*ones(1,long);     %燃速(m/s)
e_t = 0*ones(1,long);     %推进药柱已烧去肉厚(m)
e_x = 0*ones(1,long);     %续航药柱已烧去肉厚(m)
s_t = s_t0*ones(1,long);     %推进药柱燃烧面实际边长(m)
s_x = s_x0*ones(1,long);     %续航药柱燃烧面实际边长(m)
m_b = rho_p*s_t0*Lp*rb_0*ones(1,long);     %燃气生成率(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %质量流率(kg/s)
m_p_int = (phi_m*p0*At / c)*dt*ones(1,long);   %质量流率积分(kg/s)
F = 0*ones(1,long);     %推力(N)
Ab = s_t0*Lp*ones(1,long);     %燃烧面积(m^2)
Ap = Ap_t0*ones(1,long);     %通气面积(m^2)
Vg = Ap*Lp +Vg_z*ones(1,long);     %自由体积(m^3)
tau = tau_0*ones(1,long);          %虚温(K)
x = x_0*ones(1,long);        %水平位移(m)
y = y_0*ones(1,long);        %竖直位移(m)
v = v_0*ones(1,long);        %火箭速度（m/s)
theta = theta_sta*pi/180*ones(1,long);        %火箭实际角度(rad)
c_d = c_d0*ones(1,long);        %弹道系数
rho_y = rho_y0*ones(1,long);              %实际空气密度(kg/m^3)
c_s = c_s0*ones(1,long);        %声速(m/s)
Ma = 0*ones(1,long);        %马赫数
C_xon = 0.157*ones(1,long);     %阻力系数(N)
PAR = ones(15,long);        %参数记录

%循环变量
%计算循环变量
i = 1;
j = 1;
sw = zeros(1,long);        %阶段选择
swc = zeros(2,100);        %阶段变化点
n_k = 1;
%角度循环变量
ii = 1;
%续航药柱点火循环变量
jj = 1;

%燃烧周长
%s =  n_t*pi*(d_t + 2*e_t) + n_x*pi*(d_x + 2*e_x);
%各阶段通气面积
%Ap = n_t*pi*(d_t + 2*e_t)^2 / 4 + n_x*pi*(d_x + 2*e_x)^2 / 4;
%龙格库塔计算公式
%     dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;

%数值计算
%循环前常参数计算
%压强项中参数
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;


%循环计算
for theta_0 = theta_sta:theta_ste:theta_end
    %角度扫描循环
    theta(1) = theta_0*pi/180;
    
    for t_t = t_sta:t_ste:t_end
        %续航药柱点火扫描循环 
        
        %最大值预赋值
        x_max = x_0;
        y_max = y_0;
        v_max = v_0;
        i_t = 0;
        while (y(i) > 0)
            %数据记录循环
            
            if((sw(i) == 3) || (sw(i) == 15) )
                dp = 0;
                i = i + 1;
            else
                %龙格库塔逐步计算，计算出新的压强值
                %dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;
                k1 = p_a*(Ab(i) / Vg(i))*( p(i)^n_p ) - p_b*p(i) / Vg(i);
                k2 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k1/2)^n_p ) - p_b*(p(i) + dt*k1/2) / Vg(i);
                k3 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k2/2)^n_p ) - p_b*(p(i) + dt*k2/2) / Vg(i);
                k4 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k3)^n_p )  - p_b*(p(i) + dt*k3) / Vg(i);
                dp = (k1 + 2*k2 + 2*k3 + k4)/6;
                i = i + 1;
                p(i) = p(i - 1) + dp*dt;
            end

                %判断此时处于哪个阶段
            sw(i) = 1;
            if(e_t(i - 1) <= ep_t )     %推进药柱是否烧完
                sw(i) = sw(i) + 0*2;
            else
                sw(i) = sw(i) + 1*2;
            end
            if(i*dt <= t_t )        %是否到达续航药柱点火时间
                sw(i) = sw(i) + 0*4;
            else
                sw(i) = sw(i) + 1*4;
            end
            if(e_x(i - 1) <= ep_x )     %续航药柱是否烧完
                sw(i) = sw(i) + 0*8;
            else
                sw(i) = sw(i) + 1*8;
            end

            %记录阶段改变时的循环数值
            if(sw(i) ~= sw(i - 1))
                swc(1,j) = i;
                swc(2,j) = sw(i);
                j = j + 1;
            end


            %分阶段计算
            switch sw(i)
                case{1}     %第一阶段：推进药柱未烧完，未到达续航药柱点火时间，续航药柱未燃烧
                    r(i) = alpha_r*p(i)^n_p;
                    e_t(i) = e_t(i - 1) + r(i)*dt;
                    s_t(i) = n_t*pi*(d_t + 2*e_t(i));
                    e_x(i) = e_x(i - 1);
                    s_x(i) = 0;
                    Ap(i) = n_t*pi*(d_t + 2*e_t(i))^2 / 4;
                    Ab(i) = s_t(i)*Lp;
                    Vg(i) = Ap(i)*Lp + Vg_z;
                    m_b(i) = rho_p*Ab(i)*r(i);
                    m_p(i) = phi_m*p(i)*At / c;
                    F(i) = 1.64*At*(p(i)-p0) ;      
                case{3}     %第二阶段：推进药柱已烧完，未到达续航药柱点火时间，续航药柱未燃烧
                    r(i) = 0;
                    e_t(i) = e_t(i - 1);
                    s_t(i) = 0;
                    e_x(i) = e_x(i - 1);
                    s_x(i) = 0;
                    Ap(i) = n_t*pi*(Dr_t)^2 / 4;
                    Ab(i) = 0;
                    Vg(i) = Ap(i)*Lp + Vg_z;
                    m_b(i) = 0;
                    m_p(i) = 0;
                    F(i) = 0;
                case{5}     %第三阶段：推进药柱未烧完，已到达续航药柱点火时间，续航药柱未烧完
                    r(i) = alpha_r*p(i)^n_p;
                    e_t(i) = e_t(i - 1) + r(i)*dt;
                    s_t(i) = n_t*pi*(d_t + 2*e_t(i));
                    e_x(i) = e_x(i - 1) + r(i)*dt;
                    s_x(i) = n_x*pi*(d_x + 2*e_x(i));
                    Ap(i) = n_t*pi*(d_t + 2*e_t(i))^2 / 4 + n_x*pi*(d_x + 2*e_x(i))^2 / 4;
                    Ab(i) = s_t(i)*Lp + s_x(i)*Lp;
                    Vg(i) = Ap(i)*Lp + Vg_z;
                    m_b(i) = rho_p*Ab(i)*r(i);
                    m_p(i) = phi_m*p(i)*At / c;
                    F(i) = 1.64*At*(p(i)-p0) ; 
                case{7}     %第四阶段：推进药柱已烧完，已到达续航药柱点火时间，续航药柱未烧完
                    r(i) = alpha_r*p(i)^n_p;
                    e_t(i) = e_t(i - 1);
                    s_t(i) = 0;
                    e_x(i) = e_x(i - 1) + r(i)*dt;
                    s_x(i) = n_x*pi*(d_x + 2*e_x(i));
                    Ap(i) = n_t*pi*(Dr_t)^2 / 4 + n_x*pi*(d_x + 2*e_x(i))^2 / 4;
                    Ab(i) = s_x(i)*Lp;
                    Vg(i) = Ap(i)*Lp + Vg_z;
                    m_b(i) = rho_p*Ab(i)*r(i);
                    m_p(i) = phi_m*p(i)*At / c;
                    F(i) = 1.64*At*(p(i)-p0) ; 
                case{13}     %第五阶段：推进药柱未烧完，已到达续航药柱点火时间，续航药柱已烧完
                    r(i) = alpha_r*p(i)^n_p;
                    e_t(i) = e_t(i - 1) + r(i)*dt;
                    s_t(i) = n_t*pi*(d_t + 2*e_t(i));
                    e_x(i) = e_x(i - 1);
                    s_x(i) = 0;
                    Ap(i) = n_t*pi*(d_t + 2*e_t(i))^2 / 4 + n_x*pi*(Dr_x)^2 / 4;
                    Ab(i) = s_t(i)*Lp;
                    Vg(i) = Ap(i)*Lp + Vg_z;
                    m_b(i) = rho_p*Ab(i)*r(i);
                    m_p(i) = phi_m*p(i)*At / c;
                    F(i) = 1.64*At*(p(i)-p0) ; 
                case{15}     %第六阶段：推进药柱已烧完，已到达续航药柱点火时间，续航药柱已烧完
                    r(i) = 0;
                    e_t(i) = e_t(i - 1);
                    s_t(i) = 0;
                    e_x(i) = e_x(i - 1);
                    s_x(i) = 0;
                    Ap(i) = n_t*pi*(Dr_t)^2 / 4 + n_t*pi*(Dr_x)^2 / 4;
                    Ab(i) = 0;
                    Vg(i) = Ap(i)*Lp + Vg_z;
                    m_b(i) = 0;
                    m_p(i) = 0;
                    F(i) = 0;
                otherwise     %其它阶段
                    r(i) = 0;
                    e_t(i) = e_t(i - 1);
                    s_t(i) = 0;
                    e_x(i) = e_x(i - 1);
                    s_x(i) = 0;
                    Ap(i) = n_t*pi*(Dr_t)^2 / 4 + n_t*pi*(Dr_x)^2 / 4;
                    Ab(i) = 0;
                    Vg(i) = Ap(i)*Lp + Vg_z;
                    m_b(i) = 0;
                    m_p(i) = 0;
                    F(i) = 0;
            end

            %参数修正
            %逐项计算其他参数的新值
            tau(i) = tau_0-0.006328*y(i-1);              %虚温
            c_s(i) = (1.4*R_1*tau(i))^0.5;          %声速
            Ma(i) = v(i - 1)/ c_s(i);               %马赫数
            m_p_int(i) = m_p_int(i - 1) +m_p(i)*dt;      %质量流率积分
            c_d(i)= 1000*D_p*D_p/(m_k - m_p_int(i));             %弹形系数
            %阻力系数计算
            for k = 1:NUM
                if (Ma(i) <= data(1,2 ) )
                    C_xon(i) = data(1,3 );
                    break;
                elseif (Ma(i) < data(k,2 ) )
                    delta = (Ma(i) - data(k - 1,2 ) ) / (data(k,2 ) - data(k - 1,2 ) );
                    C_xon(i) = data(k - 1,3 ) + (data(k,3 ) - data(k - 1,3 ) )*delta;
                    break;
                elseif (Ma(i) == data(k,2 ) )
                    C_xon(i) = data(k,3 );
                    break;
                else
                    C_xon(i) = 0.26;
                end
            end
            %运动轨迹计算
            rho_y(i) = rho_y0*((1-2.1904e-5*y(i-1)))^4.4;                   
            v(i) = v(i-1)+ (F(i) / (m_k - m_p_int(i) ) - c_d(i)* C_xon(i)*(pi / 8000 )*v(i-1)*v(i-1)*(rho_y(i) / rho_y0 )*rho_y(i)-g*sin(theta(i-1) ) )*dt;
            theta(i) = theta(i-1)-g*cos(theta(i-1))*dt/v(i-1);
            x(i) = x(i-1)+ v(i-1)*cos(theta(i-1))*dt;
            y(i) = y(i-1)+ v(i-1)*sin(theta(i-1))*dt;
            
            %续航药柱点火记录
            if((i*dt > t_t )&&(i_t == 0) )
                i_t = i - 1;
            end

            %最大值记录
            if(x(i) > x_max)
                x_max = x(i);
            end
            if(y(i) > y_max)
                y_max = y(i);
            end
            if(v(i) > v_max)
                v_max = v(i);
            end

                %循环结束条件
            if i >= long
                break;
            end

            if( (y(i) < y(i-1) )&&(y(i) < y_end)&&(y_end < y_max) )
                break;
            end
            %数据记录循环结束
        end
         
        %计算循环后数据处理
        t_max = dt*(i - 1);     %燃烧总时间
        if((t_max < t_t )||(y_max < y_end ) )
            break;
        end
        %数据记录
        
        PAR(1,n_k) = theta_0;      %初始射角
        PAR(2,n_k) = t_max;      %运动总时间
        PAR(3,n_k) = x_max;      %射程
        PAR(4,n_k) = y_max;      %最大高度
        PAR(5,n_k) = v_max;      %最大速度
        if(y_max > y_end)   %是否到达投放高度
            PAR(6,n_k) = 1;      %已到达投放高度
            PAR(7,n_k) = v(i);      %投放时速度
            PAR(8,n_k) = theta(i)*180/pi;      %投放时角度
            if(r(i) == 0  )     %投放时药柱是否燃烧完
                PAR(9,n_k) = 0;     %投放时药柱已燃烧完
            else
                PAR(9,n_k) = 1;     %投放时药柱未燃烧完
            end
        else
            PAR(6,n_k) = 0;      %未到达投放高度
            PAR(7,n_k) = 0;      
            PAR(8,n_k) = 0;   
            PAR(9,n_k) = 0;     
        end
        if( t_max > t_t)    %续航药柱是否点火
            PAR(10,n_k) = 1;      %续航药柱已点火
            PAR(11,n_k) = t_t;      %续航药柱点火时间
            PAR(12,n_k) = x(i_t);      %续航药柱点火时水平位置
            PAR(13,n_k) = y(i_t);      %续航药柱点火时高度
            PAR(14,n_k) = v(i_t);      %续航药柱点火时速度
            PAR(15,n_k) = theta(i_t)*180/pi;      %续航药柱点火时角度
        else
            PAR(10,n_k) = 0;      %续航药柱未点火
            PAR(11,n_k) = 0;      
            PAR(12,n_k) = 0;      
            PAR(13,n_k) = 0;
            PAR(14,n_k) = 0;
            PAR(15,n_k) = 0;    
        end
        %循环参数重置
        i = 1;
        j = 1;
        sw = zeros(1,long);        %阶段选择
        swc = zeros(2,100);        %阶段变化点
        %循环计次
        n_k = n_k + 1;
        
        %续航药柱点火扫描循环结束 
        jj = jj + 1;
    end
    %循环参数重置
    jj = 1;
    i = 1;
    j = 1;
    sw = zeros(1,long);        %阶段选择
    swc = zeros(2,100);        %阶段变化点
    %角度扫描循环结束
    ii = ii + 1;
end


%循环后处理

PAR(:,n_k:long) = [];
n_k = n_k - 1;
str = [' 循环总次数 = ',num2str(n_k)];
disp(str); 


n_PAR = 1:1:n_k;

%数据输出

figure;
hold on;
plot(n_PAR,PAR(3,:));
grid on,axis tight ;
title('射程')

%写入excel
FileName = '数据记录.xlsx';
SheetName = 'Sheet1';

Str = '序号';
xlswrite(FileName,{Str},SheetName,'A1');
xlswrite(FileName,n_PAR',SheetName,'A2');

Str = '初始射角';
xlswrite(FileName,{Str},SheetName,'B1');
Str = '运动总时间';
xlswrite(FileName,{Str},SheetName,'C1');
Str = '射程';
xlswrite(FileName,{Str},SheetName,'D1');
Str = '最大高度';
xlswrite(FileName,{Str},SheetName,'E1');
Str = '最大速度';
xlswrite(FileName,{Str},SheetName,'F1');
Str = '是否到达投放高度';
xlswrite(FileName,{Str},SheetName,'G1');
Str = '投放时速度';
xlswrite(FileName,{Str},SheetName,'H1');
Str = '投放时角度';
xlswrite(FileName,{Str},SheetName,'I1');
Str = '投放时药柱是否燃烧完';
xlswrite(FileName,{Str},SheetName,'J1');
Str = '续航药柱是否点火';
xlswrite(FileName,{Str},SheetName,'K1');
Str = '点火时间';
xlswrite(FileName,{Str},SheetName,'L1');
Str = '点火时水平位置';
xlswrite(FileName,{Str},SheetName,'M1');
Str = '点火时高度';
xlswrite(FileName,{Str},SheetName,'N1');
Str = '点火时速度';
xlswrite(FileName,{Str},SheetName,'O1');
Str = '点火时角度';
xlswrite(FileName,{Str},SheetName,'P1');

xlswrite(FileName,PAR',SheetName,'B2');

