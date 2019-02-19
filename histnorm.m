function []=histnorm(T,N,n,mu,omc,lambda,H,a)
%  --- a=0,1 ---
% a=0 is used to construct the f3 histogram estimators
% a=1 is used to construct the f4 histogram estimators
sigma=sqrt(omc); 
%% Simulating random effects 
rng(231)
for r=1:10
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
y=(1/T)*X(n+1,:)+(a*lambda*D/T)*sum(X(1:n,:));
%% Histogram construction
x=[-4:0.15:5.4]; delta=0.15; p=length(x); f=zeros(p,1); 
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
    axis([-4 6 0 0.7])
    if(H ~= 0.25)
        set(gca, 'YTickLabel', '') 
    end
    if(H== 0.25)
        ylabel('density')
    end
    xlabel('\phi')
    title(['H=' num2str(H)])
    %plot
    plot([x(s); x(s+1); x(s+1)],[f(s); f(s); f(s+1)],'--','color',[0.8,0.8,0.8],'LineWidth',2)
hold on;
end
end
%% Exact density function
z=zeros(p,1);
 for s=1:p
 z(s)=1/(sqrt(2*pi)*sigma)*exp((-0.5/omc)*(x(s)-mu)^2);
 end
plot(x,z,'blue','LineWidth',2)

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

end