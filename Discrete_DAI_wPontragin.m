clear all; clc; close all;
rand('seed',1);
%% DISCRETIZATION
% B=[0 0.1065 0.0902];A=poly([1.1 0.8]);
% Gs = tf(B,A);
% a1=0;a2=0;b0=0.1;b1=0.2; %Shah's
Bp=[0 0 1];Ap=[1 1 0];Gs=tf(Bp,Ap); %Create continuous time transfer function
Ts=0.5; Hd=c2d(Gs,Ts,'matched'); % Transform continuous system to discrete system
B = Hd.Numerator{1}; A = Hd.Denominator{1};
b0=0.1; b1=0.1; a0=0.1; a1=0.01; a2=0.01;

%% RLSwEF
Am=poly([0.2+0.2j 0.2-0.2j]);Bm=[0 0.1065 0.0902];
am0=Am(1);am1=Am(2);am2=Am(3);a0=0;
Rmat=[];
factor = 25;
% Reference
T_ref = 25; t_max = 100; time = 0:Ts:t_max; nt = length(time);
% slew stuff
Tslew = 1.5; Uc = zeros(length(nt));
syms C2 t C1 c a
y(t) =C2*exp(-t) - exp(-t)*(C1*exp(t) + (t*exp(-c))/2 + (a*t*exp(t))/2);
ydot(t) = exp(-t)*(C1*exp(t) + (t*exp(-c))/2 + (a*t*exp(t))/2) - exp(-t)*(exp(-c)/2 + C1*exp(t) + (a*exp(t))/2 + (a*t*exp(t))/2) - C2*exp(-t);
for j=1:nt
    % pos or neg
    if mod(time(j),2*T_ref)<T_ref
        pn = 1;
    else
        pn =-1;
    end
    % slew
    if mod(time(j),T_ref)<Tslew

        if time(j)>=0 && time(j)<Tslew
            eqns = [y(0) == 0,ydot(0) == 0 ,y(Tslew) == 1,ydot(0) == 0];
            S = solve(eqns,[C2 C1 c a]);
            Uc(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
        else if time(j)>=50 && time(j)<50 + Tslew
                t = time(j);
                eqns = [y(50) == -1,ydot(0) == 0 ,y(50 + Tslew) == 1,ydot(0) == 0];
                S = solve(eqns,[C2 C1 c a]);

                Uc(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
        else if time(j)>=25 && time(j)<25 + Tslew
                t = time(j);
                eqns = [y(25) == 1,ydot(0) == 0 ,y(25 + Tslew) == -1,ydot(0) == 0];
                S = solve(eqns,[C2 C1 c a]);
                Uc(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
        else if time(j)>=75 && time(j)<75 + Tslew
                t = time(j);
                eqns = [y(75) == 1,ydot(0) == 0 ,y(75 + Tslew) == -1,ydot(0) == 0];
                S = solve(eqns,[C2 C1 c a]);
                Uc(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
        else if time(j)>=100 && time(j)<100 + Tslew
                t = time(j);
                eqns = [y(100) == -1,ydot(0) == 0 ,y(100 + Tslew) == 1,ydot(0) == 0];
                S = solve(eqns,[C2 C1 c a]);
                Uc(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
        end
        end
        end
        end
        end
    else
        Uc(j)=pn;
    end
end
n=4;lambda=0.95;
nzeros=2;time=zeros(1,nzeros);Y=zeros(1,nzeros);Ym=zeros(1,nzeros);
U=ones(1,nzeros);Uc=[zeros(1,nzeros),Uc];
Noise = 0;
P=[100 0 0 0;0 100 0 0;0 0 1 0;0 0 0 1]; THETA_hat(:,1)=[-a1 -a2 b0 b1]';beta=[];
alpha = 0.5; gamma = 1.2;
for i=1:nt
    phi=[]; t=i+nzeros; time(t)=i;
    Y(t)=[-A(2) -A(3) B(2) B(3)]*[Y(t-1) Y(t-2) U(t-1) U(t-2)]';
    Ym(t)=[-Am(2) -Am(3) Bm(2) Bm(3)]*[Ym(t-1) Ym(t-2) Uc(t-1) Uc(t-2)]';
    BETA=(Am(1)+Am(2)+Am(3))/(b0+b1); beta=[beta BETA];
    %RLS implementation
    phi=[Y(t-1) Y(t-2) U(t-1) U(t-2)]'; K=P*phi*1/(lambda+phi'*P*phi); P=P-P*phi*inv(1+phi'*P*phi)*phi'*P/lambda; %RLS-EF
    error(i)=Y(t)-phi'*THETA_hat(:,i); THETA_hat(:,i+1)=THETA_hat(:,i)+K*error(i);
    a1=-THETA_hat(1,i+1);a2=-THETA_hat(2,i+1);b0=THETA_hat(3,i+1);b1=THETA_hat(4,i+1);
    Af(:,i)=[1 a1 a2]'; Bf(:,i)=[b0 b1]';
    % Determine R,S, & T for CONTROLLER
    r1=(b1/b0)+(b1^2-am1*b0*b1+am2*b0^2)*(-b1+a0*b0)/(b0*(b1^2-a1*b0*b1+a2*b0^2));
    s0=b1*(a0*am1-a2-am1*a1+a1^2+am2-a1*a0)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(am1*a2-a1*a2-a0*am2+a0*a2)/(b1^2-a1*b0*b1+a2*b0^2);
    s1=b1*(a1*a2-am1*a2+a0*am2-a0*a2)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(a2*am2-a2^2-a0*am2*a1+a0*a2*am1)/(b1^2-a1*b0*b1+a2*b0^2);
    R=[1 r1];S=[s0 s1];T=BETA*[1 a0];

    Rmat=[Rmat r1];

    %calculate control signal
    U(t)=[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';
    U(t)=1.3*[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';% Arbitrarily increased to duplicate text
end
U_rlswef = U;

%% DAI
%Create command signal, Uc based on Example 3.5 plots...square wave with 50 sec period
t_max = nt-1;
THETA_hat(:,1)=[-a1 -a2 b0 b1]';
n = length(THETA_hat);
% Sigma=1/25; Noise=Sigma*randn(nt,1);
% Noise = 0;
nzeros=2;
Y_true=zeros(1,nzeros);Ym=zeros(1,nzeros);U=zeros(1,nzeros);
P=[100 0 0 0;0 100 0 0;0 0 1 0;0 0 0 1];
lambda = 1;
eb = Y_true(1) - Uc(1);
err = 0;
kp = 2.0;
kd = 6.0;
hatvec = zeros(4,1);
for i=1:t_max+1 %Loop through the output data Y(t)
    t=i+nzeros;
    de = err-eb;
    u = kp*err + kd*de;
    U(t-1) = u;
    Y_true(t)=[Y_true(t-1) Y_true(t-2) U(t-1) U(t-2)]*[-A(2) -A(3) B(2) B(3)]';
    phid = [Y_true(t) -Y_true(t-1) Y_true(t-2) -U(t-2)];
    newest = phid\u;
    hatvec(:,i) = newest;
    eb = err;
    %disp(t);
    err = Uc(t)-Y_true(t);
end
U_DAI = U;

%% PLOT
tspan = linspace(0,100,nt);
tspan = [zeros(1,2) tspan];
figure(1); %DAI
plot(tspan(1:nt),Uc(1:nt),'k-','LineWidth',1); hold on; plot(tspan(1:nt),Y_true(2:nt+1),'b--','LineWidth',3); hold off
xlabel('Time(sec)');ylabel('Output(Y)'); title('Discrete DAI using Pontryagin T = 0.35'); leg-end('Uc','Y','fontsize',11);
set(gca,'fontsize',16); set(gca,'fontname','Palatino Linotype'); xlim([0 max(time)]); grid;
% p=plot(tspan,Uc(1:203),'-',tspan,Y,'-'); p(2).LineWidth = 2; legend('Uc','Y','fontsize',11); %DAI
axis([0 100,-1.5 1.5]);
figure(2); %RLS estimation
plot(tspan(1:nt),Uc(1:nt),'k-','LineWidth',1); hold on; plot(tspan(1:nt),Y(3:nt+2),'r--','LineWidth',3); hold off
xlabel('Time(sec)');ylabel('Output(Y)'); title('RLSwEF using Pontryagin T = 0.35'); leg-end('Uc','Y','fontsize',11);
set(gca,'fontsize',16); set(gca,'fontname','Palatino Linotype'); xlim([0 max(time)]); grid;
axis([0 100,-1.5 1.5]);
DAI_err_mean = mean(abs(Uc(1:nt)-Y_true(2:nt+1)))
DAI_err_std = std(abs(Uc(1:nt)-Y_true(2:nt+1)))
RLS_err_mean = mean(abs(Uc(1:nt)-Y(3:nt+2)))
RLS_err_std = std(abs(Uc(1:nt)-Y(3:nt+2)))
Uinput_sum_RLS = mean(abs(U_rlswef))
Uinput_sum_DAI = mean(abs(U_DAI))
