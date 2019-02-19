function []=MixGamma(T,N,n,a1,b1,lambda,H,r,a)
% ---- Target: MixedGamma density function ------
u1=gamrnd(a1,b1,[2000 1]);
xx=[0:0.05:7];
%u2=gamrnd(a2,b2,[900 1]);
PHII=u1;
%+v*u2;
rng(231)
%% Kernel density estimators (25 independent samples )
for r=1:r
D=T/n; phi=gamrnd(a1,b1,[N 1]);
%+v*gamrnd(a2,b2,[N 1]);
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
y=(1/T)*X(n+1,:)+(a*lambda*D/T)*sum(X(1:n,:));
[f,xi]=ksdensity(y,xx);

%For the diagramm of the plot 
axis([0 8 0 0.7])
if(H ~= 0.25)
     set(gca, 'YTickLabel', '') 
end
if(H== 0.25)
    ylabel('density')
end
xlabel('\phi')
plot(xi,f,'color',[0.8,0.8,0.8])
hold on;
end
%% Exact density function
[f2,x2]=ksdensity(PHII,xx);
% p=length(x2); z=zeros(p,1);
%  for s=1:p
%  z(s)=1/(sqrt(2*pi)*sigma)*exp((-0.5/omc)*(x2(s)-mu)^2);
%  end
 z=gampdf([0:0.05:max(x2)],a1,b1);
 %+ gampdf(x2,a2,b2); 
plot([0:0.05:max(x2)],z,'blue','LineWidth',2) 
%% Kernel density estimate based on exact random effects
plot(x2,f2,'color',[0.4,0.4,0.4],'LineWidth',2)

if(H==0.85)
    legend([plot(xi,f,'color',[0.8,0.8,0.8])  plot([0:0.05:max(x2)],z,'blue','LineWidth',2) plot(x2,f2,'color',[0.4,0.4,0.4],'LineWidth',2)],{'Kernel estimates','True density',sprintf('Empirical random\n effects density')},'Location','northeast')
end
% xi
% x2
length(x2);
length(xi);


end