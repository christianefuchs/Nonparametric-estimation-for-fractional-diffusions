function []=PhiHat(T,N,n,mu,omc,lambda,H,s)
sigma=sqrt(omc); 
% ---- Target: Normal density function ------
PHII=sigma*randn(900,1)+mu;
%% Kernel density estimators (25 independent samples )
rng(231)
for r=1:s
D=T/n; phi=sigma*randn(1,N)+mu;
W=zeros(n,N); DW=zeros(n,N); X=zeros(n+1,N);
for i=1:N
    [w,t]=fbm1d(H,n,T);
    W(:,i)=w(2:n+1);
end
DW(1,:)=W(1,:);
for  j=2:n
    DW(j,:)=W(j,:)-W(j-1,:);
end
for i=1:N
    for j=2:n+1
        X(j,i)=X(j-1,i)+(-lambda*X(j-1,i)+phi(i))*D+DW(j-1,i);
    end
end
y=(1/T)*X(n+1,:);
[f,xi]=ksdensity(y);

%For the diagramm of the plot 
axis([-4 6 0 0.7])
plot(xi,f,'color',[0.8,0.8,0.8])
title(['H=' num2str(H)])
if(H ~= 0.25)
     set(gca, 'YTickLabel', '') 
end
if(H== 0.25)
    ylabel('density')
end
xlabel('\phi')
hold on;
end

%% Exact density function
[f2,x2]=ksdensity(PHII);
%x11=[mu-4*omc:0.1:mu+4*omc];
 p=length(x2); z=zeros(p,1);
 for s=1:p
 z(s)=1/(sqrt(2*pi)*sigma)*exp((-0.5/omc)*(x2(s)-mu)^2);
 end
plot(x2,z,'blue','LineWidth',2)

%% Kernel density estimate based on exact random effects
plot(x2,f2,'color',[0.4,0.4,0.4],'LineWidth',2)
end