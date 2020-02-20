function [Kp, Ki, K_filt_sos_num, K_filt_sos_den] = separate_I_filter_z(K) 

    %Separate controller into PI and Filter
    %The numerator of the PI is contained in the numerator of K and has to be extracted

    Ts = K.Ts;
    K_filt = K/dcgain(K);

    %The PI is scaled to obtain the initial controller
    K_PI = tf(dcgain(K), [1 -1], Ts) ;

    % Extract Kp and Ki
    K_PI =  pid(K_PI);
    Kp = K_PI.Kp;
    Ki = K_PI.Ki*Ts;

    % Split filter into SOS
    K_filt_sos = tf2sos(K_filt.num{1},K_filt.den{1},'up','inf');
    %tf2sos(K_filt.num{1},K_filt.den{1});

    K_filt_sos_den = [K_filt_sos(:,5)'; K_filt_sos(:,6)'];
    K_filt_sos_den = K_filt_sos_den(:); 

    K_filt_sos_num = K_filt_sos(:,1:3)';
    K_filt_sos_num = K_filt_sos_num(:); 
    
    if max(abs(K_filt_sos_num)) > 60
        largestCoefs = max(abs(K_filt_sos(:,1:3))');
        %largestCoefs
        K_filt_sos_num = K_filt_sos(:,1:3)';
        KsosnumRescaled = (K_filt_sos_num./(largestCoefs)).*prod(largestCoefs)^(1/length(largestCoefs));
        %KsosnumRescaled
        K_filt_sos(:,1:3) = KsosnumRescaled';
        K_filt_sos_num = K_filt_sos(:,1:3)';
        
        K_filt_sos_num = K_filt_sos_num(:); 
    end

end