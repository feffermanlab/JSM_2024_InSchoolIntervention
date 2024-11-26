clear
%Transmission rates
%Shared global common resources
betac = 0.05;
%Shared local common resources
betab = 0.1;
%Within classrooms
betaa = 0.3;

%Adjacency matrix for shared local resources
Mat=[0,1,1,0,0,0,0,0,0,0;
    1,0,1,0,0,0,0,0,0,0;
    1,1,0,0,0,0,0,0,0,0;
    0,0,0,0,1,1,1,0,0,0;
    0,0,0,1,0,1,1,0,0,0;
    0,0,0,1,1,0,1,0,0,0;
    0,0,0,1,1,1,0,0,0,0;
    0,0,0,0,0,0,0,0,1,1;
    0,0,0,0,0,0,0,1,0,1;
    0,0,0,0,0,0,0,1,1,0];
N= size(Mat,1);

%classroom sizes
sizes=[30,30,30,30,30,30,30,30,30,30];

%Relative classroom sizes
clsrm=sizes'/mean(sizes);



A =betac* ones(N,N)+(betab-betac)*Mat+(betaa-betac)*eye(N);

%Proportion of students coming to school sick
%as a constant


%time dependent phi must be defined in odefcn
phi=1;

%behavioral change from global risk 
deltac = 1;
%behavioral change from local risk
deltab = 1;
%behavioral change within classrooms
deltaa = 1;

B = deltac* ones(N,N)+(deltab-deltac)*Mat+(deltaa-deltac)*eye(N);

W=A.*B;

%recovery rate
gamma = 0.3;

%initial data
tspan = [0,100];
initiali = zeros(N,1);
initiali(1)=0.1;
y0 = [clsrm-initiali;initiali;zeros(N,1);zeros(N,1)];

[t,y]=ode78(@(t,y)odefcn(t,y,N,phi,W,gamma),tspan,y0);
size(y,1)
plot(t,y(:,N+2),t,y(:,N+1),t,y(:,2*N-1))
figure()
[phis, phires]=scanphi(0,1,200,N,W,gamma,y0,tspan);
plot(phis,phires(:,1),phis, phires(:,2),phis,phires(:,4),phis, phires(:,8))
title('Total Infection by Infected Attendence Rate')
xlabel('proportion of infected students coming to class')
ylabel('Total Infection over 100 days')
legend({'room 1','room 2','room 4', 'room 8'}, 'location','southeast')

figure()
[deltaas, deltaares]=scandeltaa(0,1.5,200,N,A,Mat,deltab,deltac,phi,gamma,y0,tspan);
plot(deltaas,deltaares(:,1),deltaas, deltaares(:,2),deltaas, deltaares(:,4),deltaas, deltaares(:,8))
title('Total Infection by within classroom behavior')
xlabel('proportion of pre-intervention behavior maintained')
ylabel('Total Infection over 100 days')
legend({'room 1','room 2','room 4', 'room 8'}, 'location','southeast')


figure()
[deltabs, deltabres]=scandeltab(0,1.5,200,N,A,Mat,deltab,deltac,phi,gamma,y0,tspan);
plot(deltabs,deltabres(:,1),deltabs, deltabres(:,2),deltabs, deltabres(:,4),deltabs, deltabres(:,8))
title('Total Infection by local behavior')
xlabel('proportion of pre-intervention behavior maintained')
ylabel('Total Infection over 100 days')
legend({'room 1','room 2','room 4', 'room 8'}, 'location','southeast')


figure()
[deltacs, deltacres]=scandeltac(0,1.5,200,N,A,Mat,deltab,deltac,phi,gamma,y0,tspan);
plot(deltacs,deltacres(:,1),deltacs, deltacres(:,2),deltacs, deltacres(:,4),deltacs, deltacres(:,8))
title('Total Infection by global behavior')
xlabel('proportion of pre-intervention behavior maintained')
ylabel('Total Infection over 100 days')
legend({'room 1','room 2','room 4', 'room 8'}, 'location','southeast')

figure()
res = scanphisense(0,1,100,N,A,Mat,gamma, y0,tspan);
plot(res(:,1),res(:,2),res(:,1),res(:,3),res(:,1),res(:,4))
title('sensitivity of total infection to interventions by $\phi$', 'interpreter','latex')
xlabel('$\phi$','interpreter','latex')
ylabel('sensitivity','interpreter','latex')
legend('$\frac{\partial J}{\partial \delta_a}$','$\frac{\partial J}{\partial \delta_b}$','$\frac{\partial J}{\partial \delta_c}$' ...
    ,'interpreter','latex','Location','southeast')


function res = scanphisense(minphi,maxphi,phires,N,A, Mat,gamma,y0,tspan)
    param = linspace(minphi,maxphi,phires);
    res = zeros(phires,4);
    for i = 1:phires
        [das ,resdas]= scandeltaa(0.998,1.002,5,N,A,Mat,1,1,param(i),gamma,y0,tspan);
        [dbs ,resdbs]= scandeltab(0.998,1.002,5,N,A,Mat,1,1,param(i),gamma,y0,tspan);
        [dcs ,resdcs]= scandeltac(0.998,1.002,5,N,A,Mat,1,1,param(i),gamma,y0,tspan);
        [M,I]=min(abs(das-1));
        h=(1.002-0.998)/5;
        partialdeltaa= sum((-resdas(I-2,:)+8*resdas(I-1,:)-8*resdas(I+1,:)+resdas(I+2,:))/(12*h));
        partialdeltab= sum((-resdbs(I-2,:)+8*resdbs(I-1,:)-8*resdbs(I+1,:)+resdbs(I+2,:))/(12*h));
        partialdeltac= sum((-resdcs(I-2,:)+8*resdcs(I-1,:)-8*resdcs(I+1,:)+resdcs(I+2,:))/(12*h));
        res(i,:)=[param(i),partialdeltaa,partialdeltab,partialdeltac];
    end
end



function [param, res]= scanphi(minphi, maxphi, phires, N,W,gamma, y0,tspan)
    param = linspace(minphi,maxphi,phires);
    res = zeros(phires, N);
    for i= 1:phires
        [t,y] = ode78(@(t,y)odefcn(t,y,N,param(i),W,gamma),tspan,y0);
        res(i,:)=trapz(t,y(:,N+1:2*N));
    end
end


function [param, res]= scandeltaa(mindeltaa, maxdeltaa, deltaares, N,A,Mat,deltab,deltac,phi,gamma,y0,tspan)
    param = linspace(mindeltaa,maxdeltaa,deltaares);
    res = zeros(deltaares, N);
    for i= 1:deltaares
        B = deltac* ones(N,N)+(deltab-deltac)*Mat+(param(i)-deltac)*eye(N);
        W = A.*B;
        [t,y] = ode78(@(t,y)odefcn(t,y,N,phi,W,gamma),tspan,y0);
        res(i,:)=trapz(t,y(:,N+1:2*N));
    end
end

function [param, res]= scandeltab(mindeltab, maxdeltab, deltabres, N,A,Mat,deltaa,deltac,phi,gamma,y0,tspan)
    param = linspace(mindeltab,maxdeltab,deltabres);
    res = zeros(deltabres, N);
    for i= 1:deltabres
        B = deltac* ones(N,N)+(param(i)-deltac)*Mat+(deltaa-deltac)*eye(N);
        W = A.*B;
        [t,y] = ode78(@(t,y)odefcn(t,y,N,phi,W,gamma),tspan,y0);
        res(i,:)=trapz(t,y(:,N+1:2*N));
    end
end

function [param, res]= scandeltac(mindeltac, maxdeltac, deltacres, N,A,Mat,deltaa,deltab,phi,gamma,y0,tspan)
    param = linspace(mindeltac,maxdeltac,deltacres);
    res = zeros(deltacres, N);
    for i= 1:deltacres
        B = param(i)* ones(N,N)+(deltab-param(i))*Mat+(deltaa-param(i))*eye(N);
        W = A.*B;
        [t,y] = ode78cd (@(t,y)odefcn(t,y,N,phi,W,gamma),tspan,y0);
        res(i,:)=trapz(t,y(:,N+1:2*N));
    end
end




function dydt = odefcn(t,y,N,phi,W,gamma)
    dydt = zeros(4*N,1);
    dydt(1:N) = -W*y(N+1:2*N).*y(1:N);
    dydt(N+1:2*N)=phi*W*y(N+1:2*N).*y(1:N)-gamma*y(N+1:2*N);
    dydt(2*N+1:3*N)=(1-phi)*W*y(N+1:2*N).*y(1:N)-gamma*y(N+1:2*N);
    dydt(3*N+1:4*N)=gamma*(y(N+1:2*N)+y(2*N+1:3*N));
end




