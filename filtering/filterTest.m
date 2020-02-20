
s = tf('s');
f1 = 2*pi*20e3;
f2 = 2*pi*100e3;
f3 = 2*pi*200e3;


csys1 = f1^2/(s^2+f1*0.3*s+f1^2);
csys2 = f2^2/(s^2+f2*0.2*s+f2^2);
csys3 = f3^2/(s^2+f3*0.1*s+f3^2);

csys = csys1 + csys2 + csys3;

dsys = c2d(csys, 1/1e6,'foh');

bode(dsys)
sosRaw = tf2sos(dsys.num{1},dsys.den{1});

b = sosRaw(:,1:3);
b = reshape(b',1,numel(b));

a = sosRaw(:,4:end);
a = reshape(a',1,numel(a));

dlmwrite('filterTest.coeffs', [a;b], ',')