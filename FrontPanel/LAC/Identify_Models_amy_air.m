% close all
clear all

data = load('amy_tube_50kHz_8191.txt','-ascii');
Ts=1/20e3;
Nprbs=12;
% data = load('amy_upper_400k_prbs11');
% Ts=1/400e3;
% Nprbs=11;

udat=data(:,1);
ydat=data(:,2); %blue position

N=2^Nprbs-1;
ini=1e3;
n=95;
u=udat(1+ini:ini+n*N);
y=ydat(1+ini:ini+n*N);

u_valid=udat((ini+5*N+1):(ini+5*N+5*N));
y_valid=ydat((ini+5*N+1):(ini+5*N+5*N));

data = detrend(iddata(y,u,Ts));
valid = detrend(iddata(y_valid,u_valid,Ts));
% plot(dat)
% DATA = iddata([],S,'Ts',Ts,'Period',per)


%% Frequency Response
U1=fft(u);
Y1=fft(y);
g=Y1(1:n:end)./U1(1:n:end);
w=0:2*pi/(Ts*N):(N-1)*2*pi/(Ts*N);
G_frd=frd(g(1:(N-1)/2),w(1:(N-1)/2));

% u_etfe = iddata([],u,'Ts',Ts,'Period',N);
% y_etfe = iddata(y,[],'Ts',Ts);
% G_etfe = etfe([y_etfe,u_etfe]);

G_spa = spa(data,700,w(1:(N-1)/2));
G_spa.ResponseData = G_spa.ResponseData;
figure;
bodemag(G_frd, G_spa)

%% ------------ Phase 2: Estimation of the order from arx loss function ------------

for i = 1 : 15 % for different possible orders
    M = arx(data, [i, i, 1]); % building an arx model
    C(i) = M.estimationInfo.LossFcn; % memorizing the value of loss fcn
end;
figure; plot(C) %order 4-9


%% ---------- Phase 3: Estimation of the order from zero-pole cancellations ---------
% We create some preliminary ARMAX models of our system and check their
% zeros-pole map.
% In the pole-zero plots you will see that at some point increasing the
% order of the model will only add more pole-zero pairs that almost cancel
% each other out. They do not improve our model, therefore we canmake a good
% guess of a sensible order of our model by looking when these pairs first appear.
% Play around with the orders and see what we get.

ordEstimLossFcn = 7; % estimated order  < C H A N G E ! ! ! >
minOrd = 3; % treating by armax system orders from ordEstimLossFcn - minOrd
plusOrd = 1; % to ordEstimLossFcn + plusOrd < C H A N G E   B O T H ! ! ! >
% zero-pole cancelation with armax
for cnt = ordEstimLossFcn - minOrd : ordEstimLossFcn + plusOrd
    cntArray = cnt + 1 - ordEstimLossFcn + minOrd; % index in the array
    MarmaxModel = armax(data, [cnt, cnt, cnt, 1]); % delay set to 1
    NoiseMod = tf(1, MarmaxModel.c, Ts); % model of the noise
    % figure; pzmap(MarmaxModel, NoiseMod); % to find the cancelations
    figure; pzmap(MarmaxModel); % to find the cancelations
    title(strcat('Order: ', num2str(cnt))) % title contains the order
    axis([-1 1 -1 1]) % to focus onto the unit circle
end;

%% ---------------- Phase 4: Estimating delay from the armax model -----------------

Marmax = armax(data, [100, 0, 0, 0]);
% estimating delay from the leading coeffs of B and their st.deviations.
% The number of coefficients that include 0 in their standard deviation is
% the number of delays
[Marmax.b + Marmax.db; Marmax.b - Marmax.db]

%% --------------- Phase 5: Estimation of the number of parameters -----------------
% We can see that at some point the error does not decrease much with an increasing number of parameter.
% Pick the order as low as possible while still having an acceptably small error
sysDelay = 1; % estimated system delay  < C H A N G E ! ! ! >
% Loss function for all possible combination
strucInt = ordEstimLossFcn - minOrd : ordEstimLossFcn + plusOrd;
nn = struc(strucInt, strucInt, sysDelay);
V = arxstruc(data, valid, nn);
close all
selstruc(V)


%% Identification and validation data
n = 6;
ni=1e5+floor(n*N/2);
zi=iddata(y(1e5:ni),u(1e5:ni),Ts);
zv=iddata(y(ni+1:(ni+1+5*N)),u(ni+1:(ni+1+5*N)),Ts);

zi=detrend(zi);
zv=detrend(zv);

% Parametric identification
na=11;nb=11;nk=1;


% Marx=arx(zi,[na nb nk]);
% Marmax=armax(zi,[na nb na nk]);
Moe=oe(zi,[nb na nk]);
% Mbj=bj(zi,[nb na na na nk]);
% Miv=iv4(zi,[na nb nk]);
% Mss=n4sid(zi);

figure;compare(zv,Moe)%,Marx,Marmax,Mbj)%,Miv,Mss)
figure;compare(G_spa,Moe)%,Marx,Marmax,Mbj)%,Miv,Mss)

%%
figure;resid(zv,Marx); legend('Arx');
figure;resid(zv,Marmax); legend('Armax');
figure;resid(zv,Moe); legend('OE');
figure;resid(zv,Mbj); legend('BJ');

model = Moe;

