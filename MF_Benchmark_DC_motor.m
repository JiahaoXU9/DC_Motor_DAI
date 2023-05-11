clear all;clc;
rand('seed',1);
Bp=[0 0 1];Ap=[1 1 0];Gs=tf(Bp,Ap); %Create continuous time transfer function
Ts=0.35; Hd=c2d(Gs,Ts,'matched'); % Transform continuous system to discrete system
B = Hd.Numerator{1}; A = Hd.Denominator{1};
a1=0;a2=0;b0=0.1;b1=0.2;
Am=poly([0.2+0.2j 0.2-0.2j]);Bm=[0 0.1065 0.0902];
am0=Am(1);am1=Am(2);am2=Am(3);a0=0;
Rmatrix=[];
factor = 1000000;
%create square wave for reference input
maxtime=200;
Uc = zeros(1,201);
for i=1:length(Uc)
    if (mod(floor(i/20),2) == 0)
        Uc(i) = 1;
    else
        Uc(i) = 0;
    end
end
traj_Uc = zeros(1,length(Uc));
check = 1;
run_next = 0;
for i=1:length(Uc)-1
    if (check)
        traj_Uc(i) = Uc(i);
        diff = Uc(i+1)-Uc(i);
        lasti = i;
        lastval = Uc(i);
    end
    if (diff ~= 0)
        check = 0;
        if (run_next)
            traj_Uc(i) = lastval + diff/2*(1+(sin(0.2*pi*(i-lasti)-pi/2)));
        end
        run_next = 1;
        if (traj_Uc(i) == Uc(i) && (i ~= lasti))
            check = 1;
            run_next = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%RECURSIVE LEAST SQUARES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uc = Uc(1:200);
n=4;lambda=1.0;
nzeros=5;time=zeros(1,nzeros);Y=zeros(1,nzeros);Ym=zeros(1,nzeros);
U=ones(1,nzeros);Uc=[ones(1,nzeros),Uc];
Noise = 1/factor*randn(1,maxtime+nzeros);
P=[100 0 0 0;0 100 0 0;0 0 1 0;0 0 0 1]; THETA_hat(:,1)=[-a1 -a2 b0 b1]';beta=[];
alpha = 0.5; gamma = 1.2;
for i=1:maxtime;
    phi=[]; t=i+nzeros; time(t)=i;
    Y(t)=[-A(2) -A(3) B(2) B(3)]*[Y(t-1) Y(t-2) U(t-1) U(t-2)]' + Noise(t-1) + Noise(t-2);
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
    Rmatrix=[Rmatrix r1];
    %calculate control signal
    U(t)=[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';
    U(t)=1.3*[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';% Arbitrarily increased to duplicate text
end
%store values of control, output, and theta for this estimation method
plotu = [U];
ploty = [Y];
plottheta = [THETA_hat];
%%%%%%%%%%%%%%%%%%%%%%%%END OF RECURSIVE LEAST SQUARES%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%AUTOREGRESSIVE MOVING AVER-AGE%%%%%%%%%%%%%%%%%%%%%
H=tf(B,A,0.5);
a1=0;a2=0;b0=0.01;b1=0.2;
Am=poly([0.2+0.2j 0.2-0.2j]);Bm=[0 0.1065 0.0902];
am0=Am(1);am1=Am(2);am2=Am(3);a0=0;
Rmatrix=[];
maxtime=200;
n=4;lambda=1;
nzeros=5;time=zeros(1,nzeros);Y=zeros(1,nzeros);Ym=zeros(1,nzeros);
U=ones(1,nzeros);
THETA_hat = zeros(4,maxtime);
THETA_hat(:,1)=[-a1 -a2 b0 b1]';beta=[];
Noise = 1/factor*randn(1,maxtime+nzeros);
epsilon=[zeros(1,nzeros+maxtime)];
n = 8;
P=10000*eye(n);P(1,1)=1000;P(2,2)=100;P(3,3)=100;P(4,4)=10000;P(5,5)=1000;P(6,6)=100;
theta_hat_els = zeros(n,1);
phi=[];
for i=1:maxtime
    t=i+nzeros; time(t)=i;
    Y(t)=[-A(2) -A(3) B(2) B(3)]*[Y(t-1) Y(t-2) U(t-1) U(t-2)]' + Noise(t-1) + Noise(t-2); %Create truth output
    BETA=(Am(1)+Am(2)+Am(3))/(b0+b1); beta=[beta BETA];
    phi=[phi; Y(t-1) Y(t-2) U(t-1) U(t-2)];
    Y(t)=[-A(2) -A(3) B(2) B(3)]*[Y(t-1) Y(t-2) U(t-1) U(t-2)]' + Noise(t-1) + Noise(t-2); %Create truth output
    BETA=(Am(1)+Am(2)+Am(3))/(b0+b1); beta=[beta BETA];
    if (i > 3)
        THETA_hat(:,i+1) = inv(phi'*phi)*phi'*Y(1+nzeros:t)';
    else
        THETA_hat(:,i+1) = THETA_hat(:,i);
    end
    a1=-THETA_hat(1,i+1);a2=-THETA_hat(2,i+1);b0=THETA_hat(3,i+1);b1=THETA_hat(4,i+1);% Up-date A & B coefficients;
    Af(:,i)=[1 a1 a2]'; Bf(:,i)=[b0 b1]'; % Store final A and B for comparison with real A&B to gener-ate epsilon errors
    % Determine R,S, & T for CONTROLLER
    r1=(b1/b0)+(b1^2-am1*b0*b1+am2*b0^2)*(-b1+a0*b0)/(b0*(b1^2-a1*b0*b1+a2*b0^2));
    s0=b1*(a0*am1-a2-am1*a1+a1^2+am2-a1*a0)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(am1*a2-a1*a2-a0*am2+a0*a2)/(b1^2-a1*b0*b1+a2*b0^2);
    s1=b1*(a1*a2-am1*a2+a0*am2-a0*a2)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(a2*am2-a2^2-a0*am2*a1+a0*a2*am1)/(b1^2-a1*b0*b1+a2*b0^2);
    R=[1 r1];S=[s0 s1];T=BETA*[1 a0];
    Rmatrix=[Rmatrix r1];
    %calculate control signal
    U(t)=[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';
    U(t)=1.3*[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';% Arbitrarily increased to duplicate text
end
plotu = [plotu; U];
ploty = [ploty; Y];
plottheta = [plottheta; THETA_hat];
%%%%%%%%%%%%%%%%%%%END OF AUTOREGRESSIVE MOVING AVER-AGE%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%EXTENDED LEAST SQUARES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=tf(B,A,0.5);
a1=0;a2=0;b0=0.01;b1=0.2;
Am=poly([0.2+0.2j 0.2-0.2j]);Bm=[0 0.1065 0.0902];
am0=Am(1);am1=Am(2);am2=Am(3);a0=0;
Rmatrix=[];
maxtime=200;
n=4;lambda=1;
nzeros=5;time=zeros(1,nzeros);Y=zeros(1,nzeros);Ym=zeros(1,nzeros);
U=ones(1,nzeros);
THETA_hat(:,1)=[-a1 -a2 b0 b1]';beta=[];% Initialize P(to), THETA_hat(to) & Beta
Noise = 1/factor*randn(1,maxtime+nzeros);
epsilon=[ones(1,nzeros+maxtime)];
n = 8;
P=10000*eye(n);P(1,1)=1000;P(2,2)=100;P(3,3)=100;P(4,4)=10000;P(5,5)=1000;P(6,6)=100;
theta_hat_els = zeros(n,1);
for i=1:maxtime;
    phi=[]; t=i+nzeros; time(t)=i;
    Y(t)=[-A(2) -A(3) B(2) B(3)]*[Y(t-1) Y(t-2) U(t-1) U(t-2)]' + Noise(t-1) + Noise(t-2); %Create truth output
    Ym(t)=[-Am(2) -Am(3) Bm(2) Bm(3)]*[Ym(t-1) Ym(t-2) Uc(t-1) Uc(t-2)]';
    BETA=(Am(1)+Am(2)+Am(3))/(b0+b1); beta=[beta BETA];
    k=i+nzeros;
    phi=[Y(t-1) Y(t-2) U(t-1) U(t-2) epsilon(t) epsilon(t-1) epsilon(t-2) epsilon(k-3)]';
    K=P*phi*1/(1+phi'*P*phi);
    P=P-P*phi*pinv(1+phi'*P*phi)*phi'*P;
    epsilon(t)=Y(t)-phi'*theta_hat_els(:,i);
    theta_hat_els(:,i+1)=theta_hat_els(:,i)+K*epsilon(t);
    THETA_hat(:,i+1) = theta_hat_els(1:4,i+1);
    a1=-THETA_hat(1,i+1);a2=-THETA_hat(2,i+1);b0=THETA_hat(3,i+1);b1=THETA_hat(4,i+1);% Up-date A & B coefficients;
    Af(:,i)=[1 a1 a2]'; Bf(:,i)=[b0 b1]'; % Store final A and B for comparison with real A&B to gener-ate epsilon errors
    r1=(b1/b0)+(b1^2-am1*b0*b1+am2*b0^2)*(-b1+a0*b0)/(b0*(b1^2-a1*b0*b1+a2*b0^2));
    s0=b1*(a0*am1-a2-am1*a1+a1^2+am2-a1*a0)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(am1*a2-a1*a2-a0*am2+a0*a2)/(b1^2-a1*b0*b1+a2*b0^2);
    s1=b1*(a1*a2-am1*a2+a0*am2-a0*a2)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(a2*am2-a2^2-a0*am2*a1+a0*a2*am1)/(b1^2-a1*b0*b1+a2*b0^2);
    R=[1 r1];S=[s0 s1];T=BETA*[1 a0];
    Rmatrix=[Rmatrix r1];
    %calculate control signal
    U(t)=[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';
    U(t)=1.3*[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';% Arbitrarily increased to duplicate text
end
plotu = [plotu; U];
ploty = [ploty; Y];
plottheta = [plottheta; THETA_hat];
%%%%%%%%%%%%%%%%%%%%%%%%%%END OF EXTENDED LEAST SQUARES%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%DETERMINISTIC AI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=tf(B,A,0.5); %Convert Plant [num] and [den] to discrete transfer function
Rmatrix=[];
%Create command signal, Uc based on Example 3.5 plots...square wave with 50 sec period
n=4;lambda=1; % number of parameters to estimate and exponential forgetting Factor
nzeros=5;time=zeros(1,nzeros);Y=zeros(1,nzeros);Ym=zeros(1,nzeros);%Initialize ouput vectors
U=ones(1,nzeros);
Noise = 1/25*randn(1,maxtime+nzeros);
epsilon=[zeros(1,nzeros+maxtime)];
n = 4;
phi_awr = [];
ustar = [];
hatvec = [];
t=[0:200];
hvy_m = [zeros(1,nzeros) traj_Uc];
eb = Y(1) - hvy_m(1);
err = 0;
kp = 2.0;
kd = 6.0;
phid = [];
ustar = [];
hatvec = zeros(4,1);
for i=1:maxtime+1; %Loop through the output data Y(t)
    t=i+nzeros; time(t)=i;
    de = err-eb;
    u = kp*err + kd*de;
    U(t-1) = u;
    Y(t)=[-A(2) -A(3) B(2) B(3)]*[Y(t-1) Y(t-2) U(t-1) U(t-2)]' + Noise(t-1) + Noise(t-2);
    phid = [phid; Y(t) -Y(t-1) Y(t-2) -U(t-2)];
    ustar = [ustar; u];
    newest = phid\ustar;
    hatvec(:,i) = newest;
    eb = err;
    err = hvy_m(t)-Y(t);
end
THETA_hat = [hatvec(2,:)./hatvec(1,:); hatvec(3,:)./hatvec(1,:); ones(1,201)./hatvec(1,:); hat-vec(4,:)./hatvec(1,:)];
plotu = [plotu; U];
ploty = [ploty; Y(1:205)];
plottheta = [plottheta; THETA_hat];
%%%%%%%%%%%%%%%%%%%%%%%%%END OF DETERMINISTIC AI%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%RECURSIVE LEAST SQUARES w/exponential forget-ting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uc = Uc(1:200);
n=4;lambda=0.99;
nzeros=5;time=zeros(1,nzeros);Y=zeros(1,nzeros);Ym=zeros(1,nzeros);
U=ones(1,nzeros);Uc=[ones(1,nzeros),Uc];
Noise = 1/factor*randn(1,maxtime+nzeros);
P=[100 0 0 0;0 100 0 0;0 0 1 0;0 0 0 1]; THETA_hat(:,1)=[-a1 -a2 b0 b1]';beta=[];
alpha = 0.5; gamma = 1.2;
for i=1:maxtime;
    phi=[]; t=i+nzeros; time(t)=i;
    Y(t)=[-A(2) -A(3) B(2) B(3)]*[Y(t-1) Y(t-2) U(t-1) U(t-2)]' + Noise(t-1) + Noise(t-2);
    Ym(t)=[-Am(2) -Am(3) Bm(2) Bm(3)]*[Ym(t-1) Ym(t-2) Uc(t-1) Uc(t-2)]';
    BETA=(Am(1)+Am(2)+Am(3))/(b0+b1); beta=[beta BETA];
    %RLS implementation
    phi=[Y(t-1) Y(t-2) U(t-1) U(t-2)]';
    K=P*phi*1/(lambda+phi'*P*phi);
    P=P-P*phi*inv(1+phi'*P*phi)*phi'*P/lambda; %RLS-EF
    error(i)=Y(t)-phi'*THETA_hat(:,i);
    THETA_hat(:,i+1)=THETA_hat(:,i)+K*error(i);
    a1=-THETA_hat(1,i+1);a2=-THETA_hat(2,i+1);b0=THETA_hat(3,i+1);b1=THETA_hat(4,i+1);
    Af(:,i)=[1 a1 a2]'; Bf(:,i)=[b0 b1]';
    % Determine R,S, & T for CONTROLLER
    r1=(b1/b0)+(b1^2-am1*b0*b1+am2*b0^2)*(-b1+a0*b0)/(b0*(b1^2-a1*b0*b1+a2*b0^2));
    s0=b1*(a0*am1-a2-am1*a1+a1^2+am2-a1*a0)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(am1*a2-a1*a2-a0*am2+a0*a2)/(b1^2-a1*b0*b1+a2*b0^2);
    s1=b1*(a1*a2-am1*a2+a0*am2-a0*a2)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(a2*am2-a2^2-a0*am2*a1+a0*a2*am1)/(b1^2-a1*b0*b1+a2*b0^2);
    R=[1 r1];S=[s0 s1];T=BETA*[1 a0];
    Rmatrix=[Rmatrix r1];
    %calculate control signal
    U(t)=[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';
    U(t)=1.3*[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';% Arbitrarily increased to duplicate text
end
%store values of control, output, and theta for this estimation method
plotu = [plotu; U];
ploty = [ploty; Y(1:205)];
plottheta = [plottheta; THETA_hat];
%%%%%%%%%%%%%%%%%%%%%%%%END OF RECURSIVE LEAST SQUARES w/exponential forgetting%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%EXTENDED LEAST SQUARES With Posterior Residuals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=tf(B,A,0.5);
a1=0;a2=0;b0=0.01;b1=0.2;
Am=poly([0.2+0.2j 0.2-0.2j]);Bm=[0 0.1065 0.0902];
am0=Am(1);am1=Am(2);am2=Am(3);a0=0;
Rmatrix=[];
maxtime=200;
n=4;lambda=1;
nzeros=5;time=zeros(1,nzeros);Y=zeros(1,nzeros);Ym=zeros(1,nzeros);
U=ones(1,nzeros);
THETA_hat(:,1)=[-a1 -a2 b0 b1]';beta=[];% Initialize P(to), THETA_hat(to) & Beta
Noise = 1/factor*randn(1,maxtime+nzeros);
epsilon=[zeros(1,nzeros+maxtime)];
n = 8;
P=10000*eye(n);P(1,1)=1000;P(2,2)=100;P(3,3)=100;P(4,4)=10000;P(5,5)=1000;P(6,6)=100;
theta_hat_els = zeros(n,1);
for i=1:maxtime;
    phi=[]; t=i+nzeros; time(t)=i;
    Y(t)=[-A(2) -A(3) B(2) B(3)]*[Y(t-1) Y(t-2) U(t-1) U(t-2)]' + Noise(t-1) + Noise(t-2); %Create truth output
    Ym(t)=[-Am(2) -Am(3) Bm(2) Bm(3)]*[Ym(t-1) Ym(t-2) Uc(t-1) Uc(t-2)]';
    BETA=(Am(1)+Am(2)+Am(3))/(b0+b1); beta=[beta BETA];
    k=i+nzeros;
    phi=[Y(t-1) Y(t-2) U(t-1) U(t-2) epsilon(t) epsilon(t-1) epsilon(t-2) epsilon(k-3)]';
    K=P*phi*1/(1+phi'*P*phi);
    P=P-P*phi*pinv(1+phi'*P*phi)*phi'*P;
    error(i)=Y(k)-phi'*theta_hat_els(:,i);
    theta_hat_els(:,i+1)=theta_hat_els(:,i)+K*error(i);
    epsilon(k)=Y(k)-phi'*theta_hat_els(:,i+1); %Form Posterior Residual
    THETA_hat(:,i+1) = theta_hat_els(1:4,i+1);
    a1=-THETA_hat(1,i+1);a2=-THETA_hat(2,i+1);b0=THETA_hat(3,i+1);b1=THETA_hat(4,i+1);% Up-date A & B coefficients;
    Af(:,i)=[1 a1 a2]'; Bf(:,i)=[b0 b1]'; % Store final A and B for comparison with real A&B to gener-ate epsilon errors
    r1=(b1/b0)+(b1^2-am1*b0*b1+am2*b0^2)*(-b1+a0*b0)/(b0*(b1^2-a1*b0*b1+a2*b0^2));
    s0=b1*(a0*am1-a2-am1*a1+a1^2+am2-a1*a0)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(am1*a2-a1*a2-a0*am2+a0*a2)/(b1^2-a1*b0*b1+a2*b0^2);
    s1=b1*(a1*a2-am1*a2+a0*am2-a0*a2)/(b1^2-a1*b0*b1+a2*b0^2)+b0*(a2*am2-a2^2-a0*am2*a1+a0*a2*am1)/(b1^2-a1*b0*b1+a2*b0^2);
    R=[1 r1];S=[s0 s1];T=BETA*[1 a0];
    Rmatrix=[Rmatrix r1];
    %calculate control signal
    U(t)=[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';
    U(t)=1.3*[T(1) T(2) -R(2) -S(1) -S(2)]*[Uc(t) Uc(t-1) U(t-1) Y(t) Y(t-1)]';% Arbitrarily increased to duplicate text
end
plotu = [plotu; U];
ploty = [ploty; Y];
plottheta = [plottheta; THETA_hat];
%%%%%%%%%%%%%%%%%%%%%%%%%%END OF EXTENDED LEAST SQUARES With Posterior Residuals %%%%%%%%%%%%%%%%%%%

y_rls = ploty(1,:);
y_arma = ploty(2,:);
y_els = ploty(3,:);
y_dai = ploty(4,:);
y_rlswef = ploty(5,:);
y_elswpr = ploty(6,:);
plot(y_rls,'g--*','LineWidth',2)
hold on
plot(y_rlswef,'r-.','LineWidth',2)
hold on
plot(y_arma,'c','LineWidth',2)
hold on
plot(y_els,'b--o','LineWidth',2)
hold on
plot(y_elswpr,':','LineWidth',2)
hold on
plot(Uc,'k-','LineWidth',2)
legend('RLS','RLSwEF','ARMA','ELS','ELSwPR','True','fontsize',11);
set(gca,'fontname','Palatino Linotype');
grid on
axis([0 55,-1.5,2.5]);
title('Comparation of different Model Following Method for DC Motor')
xlabel('Time step (in sec)'); ylabel('Output (Y)');

mean_abs_rls = mean(abs(y_rls - Uc))
mean_abs_rlswef = mean(abs(y_rlswef - Uc))
mean_abs_arma = mean(abs(y_arma - Uc))
mean_abs_els = mean(abs(y_els - Uc))
mean_abs_elswpr = mean(abs(y_elswpr - Uc))
std_abs_rls = std((y_rls - Uc))
std_abs_rlswef = std((y_rlswef - Uc))
std_abs_arma = std(y_arma - Uc)
std_abs_els = std((y_els - Uc))
std_abs_elswpr = std((y_elswpr - Uc))
mean_rls_input = mean(abs(plotu(1,:)))
mean_rlswef_input = mean(abs(plotu(2,:)))
mean_arma_input = mean(abs(plotu(3,:)))
mean_els_input = mean(abs(plotu(4,:)))
mean_elswpr_input = mean(abs(plotu(5,:)))
