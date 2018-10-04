close all; clear variables; clc

% Lab 1: Preperation. Inverted pendulum

% Variables
Mp = 0.1270;
Lr = 0.1270;
Lp = 0.3111;
Jr = 0.0083;
Jp = 0.0012;
Dr = 0.0690;
Co = 0.1285;
g = 0.981;

% A - matrix
denom = ((4*Jr + 4*Mp*Lr^2)*Jp + Mp*Lp^2*Jr);

df1_1 = 0;
df1_2 = 0;
df1_3 = 1;
df1_4 = 0;
df2_1 = 0;
df2_2 = 0;
df2_3 = 0;
df2_4 = 1;

df3_1 = 0;
df3_2 = (1/2*(20*Mp^2*Lp^2*Lr*g))/denom;
df3_3 = (-1/2*(8*Dr*Jp) - (1/2 * 2*Mp*Lp^2*Dr))/denom;
df3_4 = 0;
df4_1 = 0;
df4_2 = (20*Jr*Mp*Lp*g + 20*Mp^2*Lr^2*Lp*g)/denom;
df4_3 = (-2*Mp*Lr*Lp*Dr)/denom;
df4_4 = 0;

A_delta = [df1_1, df1_2, df1_3, df1_4;...
    df2_1, df2_2, df2_3, df2_4;...
    df3_1, df3_2, df3_3, df3_4;...
    df4_1, df4_2, df4_3, df4_4];

% B - Matrix
df1_du = 0;
df2_du = 0;
df3_du = (-1/2 * (-8*Co*Jp) + 1/2*(2*Mp*Lp^2*Co))/denom;
df4_du = (2*Mp*Lr*Lp*Co)/denom;

B_delta = [df1_du; df2_du; df3_du; df4_du];

% C - D - Matrices
C_delta = [1 0 0 0; 0 1 0 0]; 
D_delta = [0; 0];

% All matrices:
A = [A_delta, zeros(size(A_delta)); zeros(size(A_delta)), A_delta];
B = [B_delta, zeros(size(B_delta)); zeros(size(B_delta)), B_delta];
C = [C_delta, zeros(size(C_delta)); zeros(size(C_delta)), C_delta];
D = [D_delta, zeros(size(D_delta)); zeros(size(D_delta)), D_delta];

% Transfer function
system = ss(A, B, C, D);
G = tf(system);

%% Exercise 2
show_data = 0;
samples = 100;
Mp_interval = linspace(0,2*Mp,samples);
Lp_interval = linspace(0.5*Lp,1.5*Lp,samples);
Jp_interval = linspace(0,2*Jp,samples);
Co_interval = linspace(0.9*Co,1.1*Co,samples);
Mode = 'Karims way';
switch Mode
    case 'Our way'
        for i = 1:4*samples
            Mp_temp = Mp;
            Lp_temp = Lp;
            Jp_temp = Jp;
            Co_temp = Co;
            if i <= 100
                Mp_temp = Mp_interval(i);
            elseif i > 100 && i <= 200
                Lp_temp = Lp_interval(i-100);
            elseif i > 200 && i <= 300
                Jp_temp = Jp_interval(i-200);
            elseif i > 300 && i <= 400
                Co_temp = Co_interval(i-300);
            end
            denom = ((4*Jr + 4*Mp_temp*Lr^2)*Jp_temp + Mp_temp*Lp_temp^2*Jr);
            % A - Matrix
            df3_2 = (1/2*(20*Mp_temp^2*Lp_temp^2*Lr*g))/denom;
            df3_3 = (-1/2*(8*Dr*Jp_temp) - (1/2 * 2*Mp_temp*Lp_temp^2*Dr))/denom;
            df4_2 = (20*Jr*Mp_temp*Lp_temp*g + 20*Mp_temp^2*Lr^2*Lp_temp*g)/denom;
            df4_3 = (-2*Mp_temp*Lr*Lp_temp*Dr)/denom;
            
            A_delta_temp = [df1_1, df1_2, df1_3, df1_4;...
                df2_1, df2_2, df2_3, df2_4;...
                df3_1, df3_2, df3_3, df3_4;...
                df4_1, df4_2, df4_3, df4_4];
            
            % B - Matrix
            df1_du = 0;
            df2_du = 0;
            df3_du = (-1/2 * (-8*Co_temp*Jp_temp) + 1/2*(2*Mp_temp*Lp_temp^2*Co_temp))/denom;
            df4_du = (2*Mp_temp*Lr*Lp_temp*Co_temp)/denom;
            
            B_delta_temp = [df1_du; df2_du; df3_du; df4_du];
            
            % All matrices:
            A_temp = [A_delta_temp, zeros(size(A_delta_temp)); zeros(size(A_delta_temp)), A_delta_temp];
            B_temp = [B_delta_temp, zeros(size(B_delta_temp)); zeros(size(B_delta_temp)), B_delta_temp];
            
            % Transfer function
            if i <= 100
                G_tot_Mp(:,:,i) = tf(ss(A_temp, B_temp, C, D));
            elseif i > 100 && i <= 200
                G_tot_Lp(:,:,i) = tf(ss(A_temp, B_temp, C, D));
            elseif i > 200 && i <= 300
                G_tot_Jp(:,:,i) = tf(ss(A_temp, B_temp, C, D));
            elseif i > 300 && i <= 400
                G_tot_Co(:,:,i) = tf(ss(A_temp, B_temp, C, D));
            end
        end
    case 'Karims way'
        Mp_temp = ureal('Mp_temp',Mp,'Percentage',100);
        Lp_temp = ureal('Lp_temp',Lp,'Percentage',50);
        Jp_temp = ureal('Jp_temp',Jp,'Percentage',100);
        Co_temp = ureal('Co_temp',Co,'Percentage',10);
        denom = ((4*Jr + 4*Mp_temp*Lr^2)*Jp_temp + Mp_temp*Lp_temp^2*Jr);
        % A - Matrix
        
        df3_2 = (1/2*(20*Mp_temp^2*Lp_temp^2*Lr*g))/denom;
        df3_3 = (-1/2*(8*Dr*Jp_temp) - (1/2 * 2*Mp_temp*Lp_temp^2*Dr))/denom;
        df4_2 = (20*Jr*Mp_temp*Lp_temp*g + 20*Mp_temp^2*Lr^2*Lp_temp*g)/denom;
        df4_3 = (-2*Mp_temp*Lr*Lp_temp*Dr)/denom;

        A_delta_temp = [df1_1, df1_2, df1_3, df1_4;...
            df2_1, df2_2, df2_3, df2_4;...
            df3_1, df3_2, df3_3, df3_4;...
            df4_1, df4_2, df4_3, df4_4];

        % B - Matrix
        df1_du = 0;
        df2_du = 0;
        df3_du = (-1/2 * (-8*Co_temp*Jp_temp) + 1/2*(2*Mp_temp*Lp_temp^2*Co_temp))/denom;
        df4_du = (2*Mp_temp*Lr*Lp_temp*Co_temp)/denom;

        B_delta_temp = [df1_du; df2_du; df3_du; df4_du];

        % All matrices:
        A_temp = [A_delta_temp, zeros(size(A_delta_temp)); zeros(size(A_delta_temp)), A_delta_temp];
        B_temp = [B_delta_temp, zeros(size(B_delta_temp)); zeros(size(B_delta_temp)), B_delta_temp];
        system_var = ss(A_temp, B_temp, C, D);
        system_sample = usample(system_var,100);
end
if show_data == 1
    figure(1)
    subplot(2,2,1)
    sigmaplot(G_tot_Mp)
    title('Mp')
    hold all
    subplot(2,2,2)
    sigmaplot(G_tot_Lp)
    title('Lp')
    subplot(2,2,3)
    sigmaplot(G_tot_Jp)
    title('Jp')
    subplot(2,2,4)
    sigmaplot(G_tot_Co)
    title('Co')
end
% [usys, info] = ucover(system_var, system);
s = tf('s');
wiM = 1;
WiM = [wiM 0.1*wiM;0.1*wiM wiM];
Wu= diag([1/100,1/100]);
Wn = diag([1/((5/180)*pi),1/((5/180)*pi),1/((5/180)*pi),1/((5/180)*pi)]);
for i = 1:4    
    Wid(i,i) = 4/(s^2 + 4*s + 4);
end
filter1 = (10^0.2)/(s+0.01); %(10^0.2+10^(-4.25)*s)/(s+0.01);
filter2 = (10^2.7)/(s+10); %(10^2.7+10^(-6.75)*s)/(s+10);
Wp = [filter1 0 0 0;0 filter2 0 0;
    0 0 filter1 0; 0 0 0 filter2];

P = [WiM zeros(size(Wu,1),size(Wp*Wid,2)) zeros(size(Wu,1),size(Wp*Wn,2)) WiM*Wu;
    Wp*G -Wp*Wid Wp*Wn Wp*G*Wu;
    G -Wid Wn G*Wu];

%% Exercise 6

Nmeas = size(G,1);
Ncon = size(WiM*Wu,2);
[K,CL,GAM] = hinfsyn(P,Nmeas,Ncon);

