%181016
%Dashaway
%ȼ�ղ�������
%Բ��װҩ������ȼ��

%190117
%����ע��
%��ΪԲ��ҩ����ҩ����������

%201101
%��Ϊ���ܶ��ҷ�����

%201211
%�����������ⵯ��

%201229
%�޸���һЩbug������ע��

%210108
%��ε��

%210113
%ɨ��ѭ��

clear;
close all;
%��ʼֵ�趨
FileName = 'C_xon.xlsx';
[data] = xlsread(FileName,'Sheet1', 'B1:D165' );
NUM = xlsread(FileName,'Sheet1', 'A2:A2' );
%ɨ��ѭ���趨
%�Ƕ�ѭ���趨
theta_sta = 30;       %��ʼ���(deg)
theta_ste = 5;       %�������(deg)
theta_end = 75;       %�������(deg)
%���ѭ���趨
t_sta = 0;      %����ҩ������ʼʱ��(s)
t_ste = 0.5;      %����ҩ����𲽽�ʱ��(s)
t_end = 30;      %����ҩ��������ʱ��(s)

%����

%���㸨������
long = 1000000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 1e-4;      %����ֵ
%ȼ���Ҳ���
Dr_t = 0.03;       %�ƽ�ҩ��ȼ�����ھ�(m)
Dr_x = 0.03;       %����ҩ��ȼ�����ھ�(m)
At = 9.62e-4;     %��ܺ����(m^2)
D_p = 0.12;      %����(m)
m_k = 7;         %��ɵ���(kg)
Vg_z = 0;       %δװҩ�ռ����(m^3)
%װҩ����
D_t = 0.03;         %�ƽ�ҩ����ʼ�⾶(m)
d_t = 0.01;       %�ƽ�ҩ����ʼ�ھ�(m)
D_x = 0.03;         %����ҩ����ʼ�⾶(m)
d_x = 0.01;       %����ҩ����ʼ�ھ�(m)
n_t = 2;        %�ƽ�ҩ������
n_x = 2;        %����ҩ������
Lp = 0.25;      %װҩ����(m)
%��֪����
p0 = 101325;    %��ʼѹǿ(Pa)
gamma = 1.2;        %���ȱ�
rho_p = 1730;       %�ܶ�(kg/m^3)
n_p = 0.302;      %ѹǿָ��
c = 1600;       %�����ٶ�(m/s)
rb_0 = 5e-3;      %?��ʼȼ��(m/s)
alpha_r = 1.7e-4;        %?ȼ��ϵ��
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?
g = 9.78030;      %�������ٶ�(m/s^2)
R_1 = 287;        %�������峣��(J/kg/K)
rho_y0 = 1.293;       %�����ܶ�(kg/m^3)
v_e = 720;       %ȼ�������ٶ�(m/s)
%�ⵯ������
v_0 = 1;       %�����ٶ�(m/s)
x_0 = 0;           %ˮƽ��ʼλ��(m)
y_0 = 1.1;           %��ֱ��ʼλ��(m)
tau_0 = 289.1;        %��ʼ����(K) 
y_end = 400;           %��ֱ����λ��(m)

%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
s_t0 = n_t*pi*d_t;      %�ƽ�ҩ��ȼ�����ʼ�ܳ�(m)
s_x0 = n_x*pi*d_x;      %����ҩ��ȼ�����ʼ�ܳ�(m)
Ap_t0 =  n_t*pi*d_t^2 / 4;        %�ƽ�ҩ����ʼͨ�����(m^2)
Vp_t0 = (n_t*pi*(D_t^2) / 4 - Ap_t0)*Lp;       %�ƽ�ҩ��ҩ�������(m^3)
Vp_x0 = (n_x*pi*(D_x^2) / 4 - n_x*pi*d_x^2 / 4)*Lp;       %����ҩ��ҩ�������(m^3)
mp = (Vp_t0 + Vp_x0)*rho_p;     %ҩ��������(kg)
c_d0 = 1000*D_p*D_p/m_k;      %��ʼ����ϵ��
c_s0 = (1.4*R_1*tau_0)^0.5;          %��ʼ����(m/s)

%Լ������

%�ƽ��ν�������
ep_t = (D_t - d_t)/2;     %�ƽ������(m)

%�����ν�������
ep_x = (D_x - d_x)/2;     %���������(m)


%����
%�����ñ������Ѹ���ֵ��
p = p0*ones(1,long);        %ʵ��ѹǿ(Pa)
r = rb_0*ones(1,long);     %ȼ��(m/s)
e_t = 0*ones(1,long);     %�ƽ�ҩ������ȥ���(m)
e_x = 0*ones(1,long);     %����ҩ������ȥ���(m)
s_t = s_t0*ones(1,long);     %�ƽ�ҩ��ȼ����ʵ�ʱ߳�(m)
s_x = s_x0*ones(1,long);     %����ҩ��ȼ����ʵ�ʱ߳�(m)
m_b = rho_p*s_t0*Lp*rb_0*ones(1,long);     %ȼ��������(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %��������(kg/s)
m_p_int = (phi_m*p0*At / c)*dt*ones(1,long);   %�������ʻ���(kg/s)
F = 0*ones(1,long);     %����(N)
Ab = s_t0*Lp*ones(1,long);     %ȼ�����(m^2)
Ap = Ap_t0*ones(1,long);     %ͨ�����(m^2)
Vg = Ap*Lp +Vg_z*ones(1,long);     %�������(m^3)
tau = tau_0*ones(1,long);          %����(K)
x = x_0*ones(1,long);        %ˮƽλ��(m)
y = y_0*ones(1,long);        %��ֱλ��(m)
v = v_0*ones(1,long);        %����ٶȣ�m/s)
theta = theta_sta*pi/180*ones(1,long);        %���ʵ�ʽǶ�(rad)
c_d = c_d0*ones(1,long);        %����ϵ��
rho_y = rho_y0*ones(1,long);              %ʵ�ʿ����ܶ�(kg/m^3)
c_s = c_s0*ones(1,long);        %����(m/s)
Ma = 0*ones(1,long);        %�����
C_xon = 0.157*ones(1,long);     %����ϵ��(N)
PAR = ones(15,long);        %������¼

%ѭ������
%����ѭ������
i = 1;
j = 1;
sw = zeros(1,long);        %�׶�ѡ��
swc = zeros(2,100);        %�׶α仯��
n_k = 1;
%�Ƕ�ѭ������
ii = 1;
%����ҩ�����ѭ������
jj = 1;

%ȼ���ܳ�
%s =  n_t*pi*(d_t + 2*e_t) + n_x*pi*(d_x + 2*e_x);
%���׶�ͨ�����
%Ap = n_t*pi*(d_t + 2*e_t)^2 / 4 + n_x*pi*(d_x + 2*e_x)^2 / 4;
%����������㹫ʽ
%     dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;

%��ֵ����
%ѭ��ǰ����������
%ѹǿ���в���
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;


%ѭ������
for theta_0 = theta_sta:theta_ste:theta_end
    %�Ƕ�ɨ��ѭ��
    theta(1) = theta_0*pi/180;
    
    for t_t = t_sta:t_ste:t_end
        %����ҩ�����ɨ��ѭ�� 
        
        %���ֵԤ��ֵ
        x_max = x_0;
        y_max = y_0;
        v_max = v_0;
        i_t = 0;
        while (y(i) > 0)
            %���ݼ�¼ѭ��
            
            if((sw(i) == 3) || (sw(i) == 15) )
                dp = 0;
                i = i + 1;
            else
                %��������𲽼��㣬������µ�ѹǿֵ
                %dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;
                k1 = p_a*(Ab(i) / Vg(i))*( p(i)^n_p ) - p_b*p(i) / Vg(i);
                k2 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k1/2)^n_p ) - p_b*(p(i) + dt*k1/2) / Vg(i);
                k3 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k2/2)^n_p ) - p_b*(p(i) + dt*k2/2) / Vg(i);
                k4 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k3)^n_p )  - p_b*(p(i) + dt*k3) / Vg(i);
                dp = (k1 + 2*k2 + 2*k3 + k4)/6;
                i = i + 1;
                p(i) = p(i - 1) + dp*dt;
            end

                %�жϴ�ʱ�����ĸ��׶�
            sw(i) = 1;
            if(e_t(i - 1) <= ep_t )     %�ƽ�ҩ���Ƿ�����
                sw(i) = sw(i) + 0*2;
            else
                sw(i) = sw(i) + 1*2;
            end
            if(i*dt <= t_t )        %�Ƿ񵽴�����ҩ�����ʱ��
                sw(i) = sw(i) + 0*4;
            else
                sw(i) = sw(i) + 1*4;
            end
            if(e_x(i - 1) <= ep_x )     %����ҩ���Ƿ�����
                sw(i) = sw(i) + 0*8;
            else
                sw(i) = sw(i) + 1*8;
            end

            %��¼�׶θı�ʱ��ѭ����ֵ
            if(sw(i) ~= sw(i - 1))
                swc(1,j) = i;
                swc(2,j) = sw(i);
                j = j + 1;
            end


            %�ֽ׶μ���
            switch sw(i)
                case{1}     %��һ�׶Σ��ƽ�ҩ��δ���꣬δ��������ҩ�����ʱ�䣬����ҩ��δȼ��
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
                case{3}     %�ڶ��׶Σ��ƽ�ҩ�������꣬δ��������ҩ�����ʱ�䣬����ҩ��δȼ��
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
                case{5}     %�����׶Σ��ƽ�ҩ��δ���꣬�ѵ�������ҩ�����ʱ�䣬����ҩ��δ����
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
                case{7}     %���Ľ׶Σ��ƽ�ҩ�������꣬�ѵ�������ҩ�����ʱ�䣬����ҩ��δ����
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
                case{13}     %����׶Σ��ƽ�ҩ��δ���꣬�ѵ�������ҩ�����ʱ�䣬����ҩ��������
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
                case{15}     %�����׶Σ��ƽ�ҩ�������꣬�ѵ�������ҩ�����ʱ�䣬����ҩ��������
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
                otherwise     %�����׶�
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

            %��������
            %�������������������ֵ
            tau(i) = tau_0-0.006328*y(i-1);              %����
            c_s(i) = (1.4*R_1*tau(i))^0.5;          %����
            Ma(i) = v(i - 1)/ c_s(i);               %�����
            m_p_int(i) = m_p_int(i - 1) +m_p(i)*dt;      %�������ʻ���
            c_d(i)= 1000*D_p*D_p/(m_k - m_p_int(i));             %����ϵ��
            %����ϵ������
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
            %�˶��켣����
            rho_y(i) = rho_y0*((1-2.1904e-5*y(i-1)))^4.4;                   
            v(i) = v(i-1)+ (F(i) / (m_k - m_p_int(i) ) - c_d(i)* C_xon(i)*(pi / 8000 )*v(i-1)*v(i-1)*(rho_y(i) / rho_y0 )*rho_y(i)-g*sin(theta(i-1) ) )*dt;
            theta(i) = theta(i-1)-g*cos(theta(i-1))*dt/v(i-1);
            x(i) = x(i-1)+ v(i-1)*cos(theta(i-1))*dt;
            y(i) = y(i-1)+ v(i-1)*sin(theta(i-1))*dt;
            
            %����ҩ������¼
            if((i*dt > t_t )&&(i_t == 0) )
                i_t = i - 1;
            end

            %���ֵ��¼
            if(x(i) > x_max)
                x_max = x(i);
            end
            if(y(i) > y_max)
                y_max = y(i);
            end
            if(v(i) > v_max)
                v_max = v(i);
            end

                %ѭ����������
            if i >= long
                break;
            end

            if( (y(i) < y(i-1) )&&(y(i) < y_end)&&(y_end < y_max) )
                break;
            end
            %���ݼ�¼ѭ������
        end
         
        %����ѭ�������ݴ���
        t_max = dt*(i - 1);     %ȼ����ʱ��
        if((t_max < t_t )||(y_max < y_end ) )
            break;
        end
        %���ݼ�¼
        
        PAR(1,n_k) = theta_0;      %��ʼ���
        PAR(2,n_k) = t_max;      %�˶���ʱ��
        PAR(3,n_k) = x_max;      %���
        PAR(4,n_k) = y_max;      %���߶�
        PAR(5,n_k) = v_max;      %����ٶ�
        if(y_max > y_end)   %�Ƿ񵽴�Ͷ�Ÿ߶�
            PAR(6,n_k) = 1;      %�ѵ���Ͷ�Ÿ߶�
            PAR(7,n_k) = v(i);      %Ͷ��ʱ�ٶ�
            PAR(8,n_k) = theta(i)*180/pi;      %Ͷ��ʱ�Ƕ�
            if(r(i) == 0  )     %Ͷ��ʱҩ���Ƿ�ȼ����
                PAR(9,n_k) = 0;     %Ͷ��ʱҩ����ȼ����
            else
                PAR(9,n_k) = 1;     %Ͷ��ʱҩ��δȼ����
            end
        else
            PAR(6,n_k) = 0;      %δ����Ͷ�Ÿ߶�
            PAR(7,n_k) = 0;      
            PAR(8,n_k) = 0;   
            PAR(9,n_k) = 0;     
        end
        if( t_max > t_t)    %����ҩ���Ƿ���
            PAR(10,n_k) = 1;      %����ҩ���ѵ��
            PAR(11,n_k) = t_t;      %����ҩ�����ʱ��
            PAR(12,n_k) = x(i_t);      %����ҩ�����ʱˮƽλ��
            PAR(13,n_k) = y(i_t);      %����ҩ�����ʱ�߶�
            PAR(14,n_k) = v(i_t);      %����ҩ�����ʱ�ٶ�
            PAR(15,n_k) = theta(i_t)*180/pi;      %����ҩ�����ʱ�Ƕ�
        else
            PAR(10,n_k) = 0;      %����ҩ��δ���
            PAR(11,n_k) = 0;      
            PAR(12,n_k) = 0;      
            PAR(13,n_k) = 0;
            PAR(14,n_k) = 0;
            PAR(15,n_k) = 0;    
        end
        %ѭ����������
        i = 1;
        j = 1;
        sw = zeros(1,long);        %�׶�ѡ��
        swc = zeros(2,100);        %�׶α仯��
        %ѭ���ƴ�
        n_k = n_k + 1;
        
        %����ҩ�����ɨ��ѭ������ 
        jj = jj + 1;
    end
    %ѭ����������
    jj = 1;
    i = 1;
    j = 1;
    sw = zeros(1,long);        %�׶�ѡ��
    swc = zeros(2,100);        %�׶α仯��
    %�Ƕ�ɨ��ѭ������
    ii = ii + 1;
end


%ѭ������

PAR(:,n_k:long) = [];
n_k = n_k - 1;
str = [' ѭ���ܴ��� = ',num2str(n_k)];
disp(str); 


n_PAR = 1:1:n_k;

%�������

figure;
hold on;
plot(n_PAR,PAR(3,:));
grid on,axis tight ;
title('���')

%д��excel
FileName = '���ݼ�¼.xlsx';
SheetName = 'Sheet1';

Str = '���';
xlswrite(FileName,{Str},SheetName,'A1');
xlswrite(FileName,n_PAR',SheetName,'A2');

Str = '��ʼ���';
xlswrite(FileName,{Str},SheetName,'B1');
Str = '�˶���ʱ��';
xlswrite(FileName,{Str},SheetName,'C1');
Str = '���';
xlswrite(FileName,{Str},SheetName,'D1');
Str = '���߶�';
xlswrite(FileName,{Str},SheetName,'E1');
Str = '����ٶ�';
xlswrite(FileName,{Str},SheetName,'F1');
Str = '�Ƿ񵽴�Ͷ�Ÿ߶�';
xlswrite(FileName,{Str},SheetName,'G1');
Str = 'Ͷ��ʱ�ٶ�';
xlswrite(FileName,{Str},SheetName,'H1');
Str = 'Ͷ��ʱ�Ƕ�';
xlswrite(FileName,{Str},SheetName,'I1');
Str = 'Ͷ��ʱҩ���Ƿ�ȼ����';
xlswrite(FileName,{Str},SheetName,'J1');
Str = '����ҩ���Ƿ���';
xlswrite(FileName,{Str},SheetName,'K1');
Str = '���ʱ��';
xlswrite(FileName,{Str},SheetName,'L1');
Str = '���ʱˮƽλ��';
xlswrite(FileName,{Str},SheetName,'M1');
Str = '���ʱ�߶�';
xlswrite(FileName,{Str},SheetName,'N1');
Str = '���ʱ�ٶ�';
xlswrite(FileName,{Str},SheetName,'O1');
Str = '���ʱ�Ƕ�';
xlswrite(FileName,{Str},SheetName,'P1');

xlswrite(FileName,PAR',SheetName,'B2');

