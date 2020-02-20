function G_freqresp = makeFreqRespFromData(data, order_prbs, n_periods, Ts, mode)

[i,j]=size(data);
if( j > i)
   data = data' ;
end
length_prbs=2^order_prbs-1;

save data_ident data

udat=data(:,1);
ydat=data(:,2);
u=udat(1:n_periods*length_prbs);
y=ydat(1:n_periods*length_prbs);

if mode==0
  data = detrend(iddata(y,u,Ts));
  w=0:2*pi/(Ts*length_prbs):(length_prbs-1)*2*pi/(Ts*length_prbs);
  G_freqresp = spa(data,round(length_prbs/4),w(1:(length_prbs-1)/2));
elseif mode==1
  w=0:2*pi/(Ts*length_prbs):(length_prbs-1)*2*pi/(Ts*length_prbs);
  U1=fft(u);
  Y1=fft(y);
  g=Y1(1:n_periods:end)./U1(1:n_periods:end);
  G_freqresp=frd(g(1:(length_prbs-1)/2),w(1:(length_prbs-1)/2));
end

% save G_freqresp G_freqresp
% 
% [mag,phase,w] = bode(G_freqresp);
% f = (w/(2*pi))';
% mag = mag2db(mag(:)');
% phase = phase(:)';
% 
% status = 'Done!';

end