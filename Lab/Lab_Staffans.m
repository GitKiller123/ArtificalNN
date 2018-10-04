

clear all
close all
clc
Mp = 0.127;
Lr = 0.127;
Lp = 0.3111;
Jr = 0.0083;
Jp = 0.0012;
Dr = 0.069;
C0 = 0.1285;
g = 0.981;


denominator = (4*Jr+4*Mp*Lr^2)*Jp + Mp*Lp^2*Jr;
df3dtheta = 0;
df3dalpha = (10*Mp^2 * Lp^2 * Lr*g)/denominator;
df3dthetadot = (-0.5*(8*Dr*Jp)-0.5*2*Mp*Lp^2*Dr)/denominator;
df3dalphadot = 0;
df3theta = 0;
df4dtheta = 0;
df4dalpha = (20*Jr*Mp*Lp*g+20*Mp^2*Lr^2 *Lp*g)/denominator;
df4dthetadot = (-2*Mp*Lr*Lp*Dr)/denominator;
df4dalphadot = 0;

Bdelta = [0;0;(4*C0*Jp+Mp*Lp^2*C0)/denominator;(2*Mp*Lr*Lp*C0)/denominator];

Adelta = [0 0 1 0;
    0 0 0 1;
    df3dtheta df3dalpha df3dthetadot df3dalphadot;
    df4dtheta df4dalpha df4dthetadot df4dalphadot];
Cdelta = [1 0 0 0;0 1 0 0];
Ddelta = [0;0];

A1 = [Adelta zeros(size(Adelta));
    zeros(size(Adelta)) Adelta];
B1 = [Bdelta zeros(size(Bdelta));
    zeros(size(Bdelta)) Bdelta];
C1 = [Cdelta zeros(size(Cdelta));
    zeros(size(Cdelta)) Cdelta];
D1 = [Ddelta zeros(size(Ddelta));
    zeros(size(Ddelta)) Ddelta];


%Uncertainty
Mp = ureal('Mp',Mp,'Percentage',100);
Lp = ureal('Lp',Lp,'Percentage',50);
Jp = ureal('Jp',Jp,'Percentage',100);
C0 = ureal('C0',C0,'Percentage',10);


denominator = (4*Jr+4*Mp*Lr^2)*Jp + Mp*Lp^2*Jr;
df3dtheta = 0;
df3dalpha = (10*Mp^2 * Lp^2 * Lr*g)/denominator;
df3dthetadot = (-0.5*(8*Dr*Jp)-0.5*2*Mp*Lp^2*Dr)/denominator;
df3dalphadot = 0;
df3theta = 0;
df4dtheta = 0;
df4dalpha = (20*Jr*Mp*Lp*g+20*Mp^2*Lr^2 *Lp*g)/denominator;
df4dthetadot = (-2*Mp*Lr*Lp*Dr)/denominator;
df4dalphadot = 0;

Bdelta = [0;0;(4*C0*Jp+Mp*Lp^2*C0)/denominator;(2*Mp*Lr*Lp*C0)/denominator];

Adelta = [0 0 1 0;
    0 0 0 1;
    df3dtheta df3dalpha df3dthetadot df3dalphadot;
    df4dtheta df4dalpha df4dthetadot df4dalphadot];
Cdelta = [1 0 0 0;0 1 0 0];
Ddelta = [0;0];

A = [Adelta zeros(size(Adelta));
    zeros(size(Adelta)) Adelta];
B = [Bdelta zeros(size(Bdelta));
    zeros(size(Bdelta)) Bdelta];
C = [Cdelta zeros(size(Cdelta));
    zeros(size(Cdelta)) Cdelta];
D = [Ddelta zeros(size(Ddelta));
    zeros(size(Ddelta)) Ddelta];

eigsys = eig(ss(Adelta,Bdelta,Cdelta,Ddelta));
transfer = tf(ss(Adelta,Bdelta,Cdelta,Ddelta));
contr = ctrb(ss(Adelta,Bdelta,Cdelta,Ddelta));
obsvr = obsv(ss(Adelta,Bdelta,Cdelta,Ddelta));

rankcon = rank(contr)
rankobs = rank(obsvr)

sample = usample(Mp,100);

[usys,info] = ucover(ss(A,B,C,D),ss(A1,B1,C1,D1))


sigmaplot((ss(Adelta,Bdelta,Cdelta,Ddelta)))









%%

close all
Mpvec = linspace(0,2*Mp,100);
Lpvec = linspace(0.5*Lp,1.5*Lp,100);
Jpvec = linspace(0,2*Jp,100);
C0vec = linspace(0.9*C0,1.1*C0,100);

Mpvec = Mp*ones(1,100);
Lpvec = Lp*ones(1,100);
Jpvec = Jp*ones(1,100);
%C0vec = C0*ones(1,100);
c = 0

for i = 1:100

    denominator = (4*Jr+4*Mpvec(i)*Lr^2)*Jpvec(i) + Mpvec(i)*Lpvec(i)^2*Jr;
    df3dtheta = 0;
    df3dalpha = (10*Mpvec(i)^2 * Lpvec(i)^2 * Lr*g)/denominator;
    df3dthetadot = (-0.5*(8*Dr*Jpvec(i))-0.5*2*Mpvec(i)*Lpvec(i)^2*Dr)/denominator;
    df3dalphadot = 0;
    df3theta = 0;
    df4dtheta = 0;
    df4dalpha = (20*Jr*Mpvec(i)*Lpvec(i)*g+20*Mpvec(i)^2*Lr^2 *Lpvec(i)*g)/denominator;
    df4dthetadot = (-2*Mpvec(i)*Lr*Lpvec(i)*Dr)/denominator;
    df4dalphadot = 0;

    Bdelta = [0;0;(4*C0vec(i)*Jpvec(i)+Mpvec(i)*Lpvec(i)^2*C0vec(i))/denominator;(2*Mpvec(i)*Lr*Lpvec(i)*C0vec(i))/denominator];

    Adelta = [0 0 1 0;
        0 0 0 1;
        df3dtheta df3dalpha df3dthetadot df3dalphadot;
        df4dtheta df4dalpha df4dthetadot df4dalphadot];
    Cdelta = [1 0 0 0;0 1 0 0];
    Ddelta = [0;0];

    A = [Adelta zeros(size(Adelta));
        zeros(size(Adelta)) Adelta];
    B = [Bdelta zeros(size(Bdelta));
        zeros(size(Bdelta)) Bdelta];
    C = [Cdelta zeros(size(Cdelta));
        zeros(size(Cdelta)) Cdelta];
    D = [Ddelta zeros(size(Ddelta));
        zeros(size(Ddelta)) Ddelta];

    transfer(:,:,i) = tf(ss(A,B,C,D));
    
   

end

%%
transfer(:,:,1)
transfer(:,:,100)

%sigmaplot(transfer(:,:,1),transfer(:,:,100))

%%

N   = 120;
Fs  = 48e3;
Fp  = 9e3;
Ap  = 0.01;
Ast = 80;

Rp  = (10^(Ap/20) - 1)/(10^(Ap/20) + 1);
Rst = 10^(-Ast/20);

NUM = firceqrip(N,Fp/(Fs/2),[Rp Rst],'passedge');
fvtool(NUM,'Fs',Fs)




    

