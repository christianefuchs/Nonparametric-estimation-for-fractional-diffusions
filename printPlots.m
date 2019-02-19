function []= printPlots()
%General parameters for the plots
lambda=3*10^(-3);
T=100;
N=1000;
n=2^8;

%density plots with H= 0.25 or 0.75 or 0.85
figure
set(gcf,'Units','centimeters','position',[8 5 20 16])
%First plot first row
%Normaldensity mu=1, omc=0.8 and s= 50
subplot('position',[0.07 0.56 0.27 0.4])
PhiHat(T,N,n,1,0.8,lambda,0.25,50)
subplot('position',[0.37 0.56 0.27 0.4])
PhiHat(T,N,n,1,0.8,lambda,0.75,50)
subplot('position',[0.67 0.56 0.27 0.4])
PhiHat(T,N,n,1,0.8,lambda,0.85,50)
%First Plot secound row
%Gammadensity with a1=2, b1=0.9, r= 50 and a= 0 
subplot('position',[0.07 0.08 0.27 0.4])
MixGamma(T,N,n,2,0.9,lambda,0.25,50,0)
subplot('position',[0.37 0.08 0.27 0.4])
MixGamma(T,N,n,2,0.9,lambda,0.75,50,0)
subplot('position',[0.67 0.08 0.27 0.4])
MixGamma(T,N,n,2,0.9,lambda,0.85,50,0)

%Saving these as an png and pdf
saveas(gcf,'density.png')
saveas(gcf,'density.pdf')


%histogram plots with H = 0.25 or 0.75 or 0.85
figure
set(gcf,'Units','centimeters','position',[8 5 20 16])
%First plot first row
%Normaldistribution mu=1, omc=0.8 and s= 0
subplot('position',[0.07 0.56 0.27 0.4])
histnorm(T,N,n,1,0.8,lambda,0.25,0)
subplot('position',[0.37 0.56 0.27 0.4])
histnorm(T,N,n,1,0.8,lambda,0.75,0)
subplot('position',[0.67 0.56 0.27 0.4])
histnorm(T,N,n,1,0.8,lambda,0.85,0)
%First Plot secound row
%Gammadensity with a1=2, b1=0.9 and a= 0 
subplot('position',[0.07 0.08 0.27 0.4])
histGamma(T,N,n,2,0.9,lambda,0.25,0)
subplot('position',[0.37 0.08 0.27 0.4])
histGamma(T,N,n,2,0.9,lambda,0.75,0)
subplot('position',[0.67 0.08 0.27 0.4])
histGamma(T,N,n,2,0.9,lambda,0.85,0)

%Saving these as an png and pdf
saveas(gcf,'histo.png')
saveas(gcf,'histo.pdf')


