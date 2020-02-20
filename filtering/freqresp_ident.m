% close all
clear all

data = load('testID.txt','-ascii');
Ts=1/40e3; %sampling time
Nprbs=13; %prbs order

udat=data(:,1);
ydat=data(:,2); 

N=2^Nprbs-1;
ini=0; % cut 1:ini samples to remove startup transients
n=5; % number of periods
u=udat(1+ini:ini+n*N);
y=ydat(1+ini:ini+n*N);

data = detrend(iddata(y,u,Ts));
% plot(dat)

%% Frequency Response
U1=fft(u);
Y1=fft(y);
g=Y1(1:n:end)./U1(1:n:end);
w=0:2*pi/(Ts*N):(N-1)*2*pi/(Ts*N);
G_frd=frd(g(1:(N-1)/2),w(1:(N-1)/2));

windowlength = 700; %change to achieve desired result
G_spa = spa(data,windowlength,w(1:(N-1)/2));
G_spa.ResponseData = G_spa.ResponseData;
figure;
bodemag(G_frd, G_spa)