function []=histGamma(T,N,n,a1,b1,lambda,H,a)
%  --- a=0,1 ---
% a=0 is used to construct the f3 histogram estimators
% a=1 is used to construct the f4 histogram estimators
% u1=gamrnd(a1,b1,[1400 1]); u2=gamrnd(a2,b2,[1400 1]);
% PHII=u*u1+v*u2;
%% Simulating random effects 
rng(231)
for r=1:10
D=T/n; phi=gamrnd(a1,b1,[N 1]);
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
%% Histogram construction
x=[0:0.15:7.4]; delta=0.15; p=length(x); f=zeros(p,1); 
for s=1:p-1
    U=0;
    for i=1:N
        if (y(i) >= x(s)) & (y(i) < x(s+1))
        U=U+1;
        end;
    end
    f(s)=U/(N*delta);
end
for s=1:p-1
    %For the diagramm of the plot 
    axis([0 8 0 0.7])
    if(H ~= 0.25)
        set(gca, 'YTickLabel', '') 
    end
    if(H== 0.25)
        ylabel('density')
    end
    xlabel('\phi')
    %plot
    plot([x(s); x(s+1); x(s+1)],[f(s); f(s); f(s+1)],'--','color',[0.8,0.8,0.8],'LineWidth',2)
hold on;
end
end
%% Exact density function
z=gampdf([0:0.05:max(x)],a1,b1); 
plot([0:0.05:max(x)],z,'blue','LineWidth',2)

%% Exact historame
for s=1:p-1
    U=0;
    for i=1:N
        if (phi(i) >= x(s)) & (phi(i) < x(s+1))
        U=U+1;
        end;
    end
    f(s)=U/(N*delta);
end
for s=1:p-1
    plot([x(s); x(s+1); x(s+1)],[f(s); f(s); f(s+1)],'-','color',[0.4,0.4,0.4],'LineWidth',2)
hold on;
end
if(H==0.85)
    legend([plot([x(s); x(s+1); x(s+1)],[f(s); f(s); f(s+1)],'--','color',[0.8,0.8,0.8],'LineWidth',2) plot([0:0.05:max(x)],z,'blue','LineWidth',2) plot([x(s); x(s+1); x(s+1)],[f(s); f(s); f(s+1)],'-','color',[0.4,0.4,0.4],'LineWidth',2)],{sprintf('Histogram \n estimates'),'True density',sprintf('Empirical random\n effects density')})
end
end