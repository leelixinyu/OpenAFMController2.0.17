function [K,status] = controldesign_hsScanner_sweep

load hsScannerSweepsSantiAmp_00000

w = dev2062.demods.sample{4}.frequency*2*pi;
G_sweep = frd(dev2062.demods.sample{4}.x + 1i*dev2062.demods.sample{4}.y, w);
G_sweep.ResponseData = G_sweep.ResponseData./abs(G_sweep.ResponseData(:,:,1));

Fs = 500e3; %Controller Sampling Frequency
Ts = 1/Fs;

s = tf('s');
z = tf('z',Ts);

%% Plant
wmin=w(1);
i = find(w < 1/Ts*pi,1,'last');
wmax= w(i); 
wobj = w(1:i);
wstab = logspace(log10(wmin),log10(1/Ts*pi),500);

G{1} = frd(G_sweep.ResponseData(:,:,1:i), G_sweep.Frequency(1:i));
G{1}.Ts = Ts;
nmod=1; %number of models
n = size(G{1},1);

%% initial controller 
Kinit = tf(1e-3);

%% Controller Structure and Initial Controller
nx = 11; ny = 10; %Controller Order

for j=1:n
    Kden = [pole(Kinit(j,j))', 0*ones(1, ny+1-length(Kinit(j,j).den{1}))] ;
    Y_c(j,j,:) = flip(poly(Kden));
    Knum = flip(Kinit(j,j).num{1} * poly(0*ones(1,nx))) ;
    X_c(j,j,:) = dcgain(Kinit(j,j))/(sum(Knum)/sum(poly(Kden)))*Knum;   
end

% Fx, Fy are multiplied element-wise!
% [b,a] = butter(2,0.95); %add lowpass filter
Fy = tf(ones(n));%*(z-1)*tf(a,1,Ts);
Fx = tf(ones(n));%*tf(b,1,Ts);

% load K_cl
% X_c = X_cl{end}; Y_c = Y_cl{end};


%% Control Performance
%Desired open-loop transfer function Ld
Ld = 1/( (1/200e3)*s+1);%50e3*2*pi/s;

%Objective Weight
clear Wobj
tmp = ones(1, length(wobj));
for i=1:length(wobj)
    Wobj(:,:,i) = tmp(i)*diag(ones(1,n));
end
Wobj = frd(Wobj,wobj);

% Sensitivity Weight
clear W1
tmp = ones(1, length(wobj));
for i=1:length(wobj)
    W1(:,:,i) = inv( tmp(i)*diag(ones(1,n)) );
end
W1 = frd(W1,wobj);

% Complementary Sensitivity Weight
clear W2
tmp = freqresp(feedback(Ld,1), wobj);
for i=1:length(wobj)
    W2(:,:,i) = inv( 1.2*tmp(i)* diag(ones(1,n)) );
end
W2 = frd(W2,wobj);

% Input Sensitivity Weight
clear W3
wc = 2.3e5*2*pi;
[b,a] = butter(2, wc/(pi/Ts));
tmp = freqresp(tf(b,a,Ts), wobj);
% tmp = tmp.*freqresp(0.1*s/(0.1*s+1),wobj); %0.02
tmp = abs(tmp(:));
tmp(end) = tmp(end-1); %last value would be 0 and cant be inverted

for i=1:length(wobj)
    W3(:,:,i) = inv( 10*tmp(i) *eye(n) );
end
W3 = frd(W3,wobj);


%% Generate Frequency Responses
for i=1:nmod
    G_w{i} = freqresp(G{i},wobj);
end

Ld_w = freqresp(Ld,wobj);
Wobj_w = freqresp(Wobj,wobj);
W1_w = freqresp(W1,wobj);
W2_w = freqresp(W2,wobj);
W3_w = freqresp(W3,wobj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Iterative Algorithm
tic
obj=0;
maxiter=2; % max number of iterations
X_cl={}; Y_cl={};

for k=1:maxiter
    
    %--------Optimization Variables-------------------------
    % New controller variables
    X = sdpvar(n,n,nx+1,'full');
    Y = sdpvar(n,n,ny+1,'diagonal');
    Y(:,:,ny+1)=eye(n);
    
    Constraints = [];   
    
    % Define objective function
    % Minimize Hinf norm between L and Ld
    Objective=0;
    %Hinf
        gamma_LS = sdpvar(n,n, nmod*length(wobj),'hermitian');
        for j=1:nmod*length(wobj)
            Objective = Objective + wobj(j)/wobj(end)*trace(gamma_LS(:,:,j)); %maximize gamma => minimize -gamma
        end

%         Constraints = [Constraints, gamma_LS > 0];        
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRIDDING CONSTRAINTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for jj=1:nmod

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Stability Constraints        
%         for j=1:length(wstab) 
% %             zj = freqresp(s,wstab(j)); 
%             zj = freqresp(z,wstab(j)); %Bilinear Transform
%             Yk=zeros(n); Yk_c = zeros(n);
%             for i=1:(ny+1)
%                 Yk = Yk + Y(:,:,i)*zj^(i-1);
%                 Yk_c = Yk_c + Y_c(:,:,i)*zj^(i-1);
%             end            
%             Yk = Yk.*evalfr(Fy,zj); Yk_c = Yk_c.*evalfr(Fy,zj); 
%             
%             %Stability Constraint (det(Y)~=0 <-> YY* > 0). Automatically guaranteed by 
%             %Loop Shaping objective function constraint
%             tmp2 = ctranspose(Yk_c)*Yk_c; % Make sure tmp = Yk_c^*Yk_c really is hermitian
%             tmp = diag(real(diag(tmp2)));
%             tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));
%             
%             Constraints = [Constraints, [ctranspose(Yk)*Yk_c+ctranspose(Yk_c)*Yk-tmp] > 0 ];
%         end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Performance Constraints    
        for j = 1:length(wobj) % Loop for constraints at each frequency point
            
            zj = freqresp(z,w(j)); %Bilinear Transform
            Xk=zeros(n); Xk_c = zeros(n);
            Yk=zeros(n); Yk_c = zeros(n);
            for i=1:(nx+1)
                Xk = Xk + X(:,:,i)*zj^(i-1);
                Xk_c = Xk_c + X_c(:,:,i)*zj^(i-1);
            end
            for i=1:(ny+1)
                Yk = Yk + Y(:,:,i)*zj^(i-1);
                Yk_c = Yk_c + Y_c(:,:,i)*zj^(i-1);
            end
            Xk = Xk.*evalfr(Fx,zj); Xk_c = Xk_c.*evalfr(Fx,zj);            
            Yk = Yk.*evalfr(Fy,zj); Yk_c = Yk_c.*evalfr(Fy,zj);            

            Gk = G_w{jj}(:,:,j); %plant evaluated at freq            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %Objective Function
            Wobjk = Wobj_w(:,:,j);    
            scalem=1;
            %Minimize ||S||_inf
%             Constraints = [Constraints, ([(ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp)./scalem(j), ctranspose(Wobjk*Yk)./sqrt(scalem(j)); ... 
%                             Wobjk*Yk./sqrt(scalem(j)), gamma_S*eye(n)] >=  0) ];    
            %Minimize ||L-Ld||_2
            Ldk = Ld_w(:,:,j);  
            tmp2 = ctranspose(Yk_c)*Yk_c;
            tmp = diag(real(diag(tmp2)));
            tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));
            scalem = abs(tmp);
                Constraints = [Constraints, [ctranspose(Yk)*Yk_c+ctranspose(Yk_c)*Yk-tmp, ctranspose(Gk*Xk-Ldk*Yk); ...
                                         Gk*Xk-Ldk*Yk, gamma_LS(:,:,j + (jj-1)*length(w))] >=0 ];
%             tmp_db(:,:,j) = [(ctranspose(Yk_c)*Yk_c+ctranspose(Yk_c)*Yk_c-tmp)/scalem, (ctranspose(Gk*Xk_c-Ldk*Yk_c))/sqrt(scalem); ...
%                                      (Gk*Xk_c-Ldk*Yk_c)/sqrt(scalem), 1];


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                     
            Pk = Yk + Gk*Xk;
            Pk_c = Yk_c + Gk*Xk_c;
            tmp2 = ctranspose(Pk_c)*Pk_c; %Make sure tmp = Pk_c^*Pk_c really is hermitian
            tmp = diag(real(diag(tmp2)));
            tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));
            
            %Sensitivity Constraint    
            W1k = W1_w(:,:,j);       
            scalem = 1;
%             Constraints = [Constraints, ([(ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp)/scalem, ctranspose(W1k*Yk)/sqrt(scalem); ... 
%                             W1k*Yk/sqrt(scalem), eye(n)] >=  0) ]; 

            %Complementary Sensitivity Constraint
            W2k = W2_w(:,:,j);
%             Constraints = [Constraints, ([(ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp)/scalem, ctranspose(W2k*Gk*Xk)/sqrt(scalem); 
%                             W2k*Gk*Xk/sqrt(scalem), eye(n)] >= 0 )];        

            %Input Sensitivity Constraint
            W3k = W3_w(:,:,j);
            Constraints = [Constraints, ([(ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp)/scalem, ctranspose(W3k*Xk)/sqrt(scalem); 
                            W3k*Xk/sqrt(scalem), eye(n)] >= 0 )];   
        end        
    end
     % Set some options for YALMIP and solver
    options = sdpsettings('verbose',2,'solver','mosek'); % mosek sedumi
    % Solve the problem
    sol = optimize(Constraints,Objective,options)
    
    % Analyze error flags
    if sol.problem == 1
     display('Hmm, something went wrong!');
     sol.info
     yalmiperror(sol.problem)
     break;
    end
    
    display(k) %to see progress...
    if(false)%abs(obj-value(Objective))/value(Objective)<0.05) % set convergence criterion
        display('Converged!')
        break; %convergence!
    else
        obj = value(Objective)        
        X_c = value(X);
        Y_c = value(Y);
        X_cl{k} = value(X);
        Y_cl{k} = value(Y);
%         save K_cl X_cl Y_cl       
    end
end
toc

%% Final Controller & Analysis
X = tf(zeros(n,n)); Yinv = tf(zeros(n,n));
for j=1:(ny+1)
    Yinv = Yinv + Y_c(:,:,j)*z^(j-1);
end
Yinv = inv(Yinv.*Fy);

for j=1:(nx+1)
    X = X + X_c(:,:,j)*z^(j-1);
end
X = X.*Fx;
K = X*Yinv;%.*[1/Fy 1; 1 1/Fy]; 

status = '.';
if sol.problem == 0
    status = 'Converged!';
elseif sol.problem == 1
    status = 'Infeasible!';
end

%%  
% Evaluate Constraints
% G_d = frd(G{1},wobj);
% G_d.Ts = Ts;
% Kfrd = frd(K, wobj);
% % load G_d_DAPI
% % Characteristic TFs
% S={};T={};U={};
% for i=1:nmod
%     S{i}=feedback(eye(n),G_d*Kfrd);
%     T{i}=feedback(G_d*Kfrd,eye(n));
%     U{i}=feedback(Kfrd,G_d);
% end
% % Plot Sensitivities
% for i=1:nmod
%     
%     figure(1)
%     W1.Ts=Ts;
%     sigma(S{i},{wmin,wmax})
%     hold on
% %     sigma(inv(W1),{wmin,wmax})
%     title('S')
%      
%     figure(2)
%     sigma(T{i},{wmin,wmax})
%     hold on
%     sigma(inv(W2),'y--')
%     sigma(feedback(Ld,1))
%     title('T')
%     
%     figure(3)
%     W3.Ts = Ts;
%     sigma(U{i},{wmin,wmax})
%     hold on
%     sigma(inv(W3),{wmin,wmax})
%     title('U')
%     
%     figure(4)
%     bodemag(G_d*Kfrd);
%     hold on
%     bodemag(Ld,wobj)
% end

end