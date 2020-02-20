

samplingFreq = 100e3;
stopFreq = 12e3;
stopFreq2 = 10e3;
Ts = 1/samplingFreq;

[b,a] = cheby2(2, 40, stopFreq/(samplingFreq/2))
[b2,a2] = butter(4, stopFreq2/(samplingFreq/2))
%%
%K=tf(b,a,Ts)*tf(b2,a2,Ts);
K=tf(b2,a2,Ts)

p = bodeoptions()
p.FreqUnits='kHz';
p.MagLowerLimMode = 'manual';
p.MagLowerLim=-50;
bodeplot(K,p)
