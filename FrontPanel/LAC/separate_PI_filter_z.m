function [Kp, Ki, K_filt_sos_num, K_filt_sos_den] = separate_PI_filter(K) 

    %Separate controller into PI and Filter
    %The numerator of the PI is contained in the numerator of K and has to be extracted

    Ts = K.Ts;
    
    %PI numerator is a first-order transfer function with a single real root
    roots_num=roots(K.num{1});
    i = find( imag(roots_num)==0,1);
    num_PI = poly(roots_num(i));

    %Filter numerator has all the other roots
    roots_num(i) = [];
    num_Kfilt = poly(roots_num);

    %The filter is scaled to have a gain of 10
    K_filt = tf(num_Kfilt, K.den{1},Ts);
    K_filt = K_filt/dcgain(K_filt);

    %The PI is scaled to obtain the initial controller
    scalePI = dcgain(tf(num_PI, 1, Ts));
    K_PI = tf(num_PI, [1 -1],Ts) * dcgain(K)/scalePI ;

    % Extract Kp and Ki
    K_PI =  pid(K_PI);
    Kp = K_PI.Kp;
    Ki = K_PI.Ki*Ts;

    % Split filter into SOS
    K_filt_sos = tf2sos(K_filt.num{1},K_filt.den{1},'up','inf');
    %tf2sos(K_filt.num{1},K_filt.den{1});

    K_filt_sos_den = [K_filt_sos(:,5)'; K_filt_sos(:,6)'];
    K_filt_sos_den = K_filt_sos_den(:); 

    K_filt_sos_num = [K_filt_sos(:,1)'; K_filt_sos(:,2)'; K_filt_sos(:,3)'];
    K_filt_sos_num = K_filt_sos_num(:); 

end