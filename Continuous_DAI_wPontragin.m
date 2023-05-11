clear all;clc;close all;
rand('seed',1);
% Enter Given Plant parameters
for k=1:2
    Bp=[0 0 1];Ap=[1 1 0];Gs=tf(Bp,Ap); %Create continuous time transfer function
    Ts=[0.5 0.35]; Hz=c2d(Gs,Ts(k),'matched'); % Transform continuous system to discrete system
    B = Hz.Numerator{1}; A = Hz.Denominator{1};
    % Initial estimates of plant parameters for undetermined system from example 3.5
    b0=0.1; b1=0.1; a0=0.1; a1=0.01; a2=0.01;
    % Reference
    T_ref = 25; t_max = 100; time = 0:Ts(k):t_max; nt = length(time);
    % slew stuff
    syms C2 t C1 c a
    Tslew = 1.5; Yd = zeros(length(nt));
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
                Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
            else if time(j)>=50 && time(j)<50 + Tslew
                    t = time(j);
                    eqns = [y(50) == -1,ydot(0) == 0 ,y(50 + Tslew) == 1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);

                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
            else if time(j)>=25 && time(j)<25 + Tslew
                    t = time(j);
                    eqns = [y(25) == 1,ydot(0) == 0 ,y(25 + Tslew) == -1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);
                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
            else if time(j)>=75 && time(j)<75 + Tslew
                    t = time(j);
                    eqns = [y(75) == 1,ydot(0) == 0 ,y(75 + Tslew) == -1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);
                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
            else if time(j)>=100 && time(j)<100 + Tslew
                    t = time(j);
                    eqns = [y(100) == -1,ydot(0) == 0 ,y(100 + Tslew) == 1,ydot(0) == 0];
                    S = solve(eqns,[C2 C1 c a]);
                    Yd(j)=S.C2.*exp(-time(j)) - exp(-time(j)).*(S.C1.*exp(time(j)) + (time(j).*exp(-S.c))/2 + (S.a.*time(j).*exp(time(j)))/2);
            end
            end
            end
            end
            end
        else
            Yd(j)=pn;
        end
    end
    THETA_hat(:,1)=[-a1 -a2 b0 b1]';
    n = length(THETA_hat);
    Sigma=1/12*0; Noise=Sigma*randn(nt,1);
    nzeros=2;Y=zeros(1,nzeros);Y_true=zeros(1,nzeros);
    Ym=zeros(1,nzeros);U=zeros(1,nzeros);Yd=[zeros(1,nzeros),Yd];
    P=[100 0 0 0;0 100 0 0;0 0 1 0;0 0 0 1];
    lambda = 1;
    for i=1:nt-1
        t=i+nzeros;
        % Update Dynamics
        Y_true(t)=[Y(t-1) Y(t-2) U(t-1) U(t-2)]*[-A(2) -A(3) B(2) B(3)]';
        Y(t)=Y_true(t)+Noise(i);
        phi=[Y(t-1) Y(t-2) U(t-1) U(t-2)]';
        K=P*phi*1/(lambda+phi'*P*phi);
        P=P-P*phi/(1+phi'*P*phi)*phi'*P/lambda;
        innov_err(i)=Y(t)-phi'*THETA_hat(:,i);
        THETA_hat(:,i+1)=THETA_hat(:,i)+K*innov_err(i);
        a1=-THETA_hat(1,i+1);a2=-THETA_hat(2,i+1);b0=THETA_hat(3,i+1);b1=THETA_hat(4,i+1);%
        THETA=[-a1 -a2 b0 b1];
        % Calculate Model control, U(t) optimally
        U(t)=[Yd(t+1) Y(t) Y(t-1) U(t-1)]*[1 a1 a2 -b0]'/b1;
    end
    Y_true(end+1)=Y_true(end);
    FS = 2;
    time = [-(nzeros-1)*Ts:Ts:0 time];
    Ts(k)
    DAI_err_mean = mean(abs(Yd-Y_true))
    DAI_err_std = std(abs(Yd-Y_true))
    U_input = mean(abs(U))
    if k==1
        figure (k)
        h1 = plot(time,Y_true,'r--','LineWidth',3);hold on;
        plot(time,Yd,'k-','LineWidth',1);
        axis([0 100,-1.5 1.5]); hold off; grid;
        legend('Uc','Y','fontsize',11); xlabel('Time(sec)');ylabel('Output(Y)');title('Continuous DAI with Pontryagin T = 0.50s'); set(gca,'fontsize',16);
        set(gca,'fontname','Palatino Linotype');
    else
        figure (k)
        h1 = plot(time,Y_true,'b--','LineWidth',3);hold on;
        plot(time,Yd,'k-','LineWidth',1);
        axis([0 100,-1.5 1.5]); hold off; grid;
        legend('Uc','Y','fontsize',11); xlabel('Time(sec)');ylabel('Output(Y)');title('Continuous DAI with Pontryagin T = 0.35s'); set(gca,'fontsize',16);
        set(gca,'fontname','Palatino Linotype');
    end
end
