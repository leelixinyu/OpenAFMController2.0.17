
function [K,Ld,status] = controldesign_z_fun_I(G_freqresp, Ts, fc, order_x, order_y, environment)

s = tf('s');
z = tf('z',Ts);

% Replace low frequencies to mitigate the effect of noise and get correct DC gain
% i = find(G_freqresp.Frequency > 0.05*1/Ts*pi,1);
% G_freqresp.ResponseData(1:i) = abs(G_freqresp.ResponseData(i));


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wobj = logspace(log10(0.01*1/Ts*pi),log10(1/Ts*pi),500);

wmin = wobj(1);
wmax = wobj(end);

G{1} = G_freqresp;
G{1}.Ts=Ts;
nmod=1; %number of models
n = size(G{1},1);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONTROLLER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if real(freqresp(G{1}, wmin)) < 0
    Kinit = tf(-1e-3);
else 
    Kinit = tf(1e-3);
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE CONTROLLER STRUCTURE, FILL IN INITIAL CONTROLLER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nx = order_x; ny = order_y; %Controller Order

for j=1:n
    Y_c(j,j,:) = [zeros(1,ny) 1];
    X_c(j,j,:) = [zeros(1,nx) dcgain(Kinit)];
end

% Fx, Fy are multiplied element-wise!
Fy = (z-1)*tf(ones(n));%*(z-1)*tf(a,1,Ts);
Fx = tf(ones(n));%*tf(b,1,Ts);

% load K_cl
% X_c = X_cl{end}; Y_c = Y_cl{end};


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE CONTROL PERFORMANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Desired open-loop transfer function Ld
Ld = fc*2*pi/s;%50e3*2*pi/s;

%Objective Weight
% clear Wobj
% tmp = ones(1, length(wobj));
% for i=1:length(wobj)
%     Wobj(:,:,i) = tmp(i)*diag(ones(1,n));
% end
% Wobj = frd(Wobj,wobj);

% Sensitivity Weight
% clear W1
% tmp = ones(1, length(wobj));
% for i=1:length(wobj)
%     W1(:,:,i) = inv( tmp(i)*diag(ones(1,n)) );
% end
% W1 = frd(W1,wobj);

% Complementary Sensitivity Weight
clear W2
wc = 1.1*2*pi*fc;
[b,a] = butter(2, wc/(pi/Ts));
tmp = freqresp(tf(b,a,Ts), wobj);
tmp = abs(tmp(:));
tmp(end) = tmp(end-1); %last value would be 0 and cant be inverted
for i=1:length(wobj)
    W2(:,:,i) = inv( 1.2*tmp(i)* diag(ones(1,n)) );
end
W2 = frd(W2,wobj);

% Input Sensitivity Weight
clear W3
wc = 1.1*fc*2*pi; %0.9*pi/Ts;
[b,a] = butter(2, wc/(pi/Ts));
tmp = freqresp( 1/abs(freqresp(G{1},wmin))*tf(b,a,Ts), wobj);
tmp = abs(tmp(:));
tmp(end) = tmp(end-1); %last value would be 0 and cant be inverted
for i=1:length(wobj)    
    if wobj(i) < fc*2*pi
        W3(:,:,i) = inv( 10*tmp(i) *eye(n) );
    else
        if environment==0 %if scanner operates in air
            W3(:,:,i) = inv( 1*tmp(i) *eye(n) ); 
        elseif environment==1 %if scanner operates in liquid
            W3(:,:,i) = inv( 1*tmp(i) *eye(n) ); 
        end
    end
end
W3 = frd(W3,wobj);


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE FREQUENCY RESPONSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:nmod
    G_w{i} = freqresp(G{i},wobj);
end
G_w{end}(end)=G_w{end}(end-1);

Ld_w = freqresp(Ld,wobj);
% Wobj_w = freqresp(Wobj,wobj);
% W1_w = freqresp(W1,wobj);
W2_w = freqresp(W2,wobj);
W3_w = freqresp(W3,wobj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATIVE ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
obj=0;
maxiter=10; % max number of iterations
X_cl={}; Y_cl={};

for k=1:maxiter
    
    %--------Optimization Variables-------------------------
    % New controller variables
    X = sdpvar(n,n,nx+1,'full');
    Y = sdpvar(n,n,ny+1,'diagonal');
    Y(:,:,ny+1)=eye(n);
    Constraints = [];
    
    % Define objective function
    Objective=0;    
    % Minimize infinity-norm of W1S
%     gamma_S = sdpvar;
%     Objective = Objective + gamma_S;
%     Constraints = [Constraints, gamma_S > 0];  

    % Minimize H2 norm between L and Ld        
    gamma_LS = sdpvar(n,n, nmod*length(wobj),'hermitian');
    for j=1:nmod*length(wobj)
        Objective = Objective + trace(gamma_LS(:,:,j));
        Constraints = [Constraints, gamma_LS(:,:,j) >= 0]; 
    end

         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRIDDING CONSTRAINTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for jj=1:nmod
        for j = 1:length(wobj) % Loop for constraints at each frequency point

            zj = freqresp(z,wobj(j)); %Bilinear Transform
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
            Yk_c_stab=Yk_c;Yk_stab=Yk;
            Xk = Xk.*evalfr(Fx,zj); Xk_c = Xk_c.*evalfr(Fx,zj);            
            Yk = Yk.*evalfr(Fy,zj); Yk_c = Yk_c.*evalfr(Fy,zj);            

            Gk = G_w{jj}(:,:,j); %plant evaluated at freq     
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FIRST STABILITY CONSTRAINT
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tmp2 = ctranspose(Yk_c_stab)*Yk_c_stab; % Make sure tmp = Yk_c^*Yk_c really is hermitian
            tmp = diag(real(diag(tmp2)));
            tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));

            Constraints = [Constraints, [ctranspose(Yk_stab)*Yk_c_stab+ctranspose(Yk_c_stab)*Yk_stab-tmp] >= 1e-10 ];   
            
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % OBJECTIVE FUNCTION
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %Minimize ||S||_inf
            %Wobjk = Wobj_w(:,:,j);    
            %Constraints = [Constraints, ([(ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp)./scalem(j), ctranspose(Wobjk*Yk)./sqrt(scalem(j)); ... 
            %                 Wobjk*Yk./sqrt(scalem(j)), gamma_S*eye(n)] >=  0) ];    

            %Minimize ||L-Ld||_2
            Ldk = Ld_w(:,:,j);  
            tmp2 = ctranspose(Yk_c)*Yk_c;
            tmp = diag(real(diag(tmp2)));
            tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));
            Constraints = [Constraints, [ctranspose(Yk)*Yk_c+ctranspose(Yk_c)*Yk-tmp, ctranspose(Gk*Xk-Ldk*Yk); ...
                                     Gk*Xk-Ldk*Yk, gamma_LS(:,:,j + (jj-1)*length(wobj))] >= 0 ];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % PERFORMANCE  CONSTRAINTS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
            Pk = Yk + Gk*Xk;
            Pk_c = Yk_c + Gk*Xk_c;
            tmp2 = ctranspose(Pk_c)*Pk_c; %Make sure tmp = Pk_c^*Pk_c really is hermitian
            tmp = diag(real(diag(tmp2)));
            tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));
            scalem = abs(tmp);
            
            % Additional stability constraint
            Constraints = [Constraints, [ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp]/scalem >= 1e-10 ];
            
            %Sensitivity Constraint    
%             W1k = W1_w(:,:,j);       
%             Constraints = [Constraints, ([(ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp)/scalem, ctranspose(W1k*Yk)/sqrt(scalem); ... 
%                             W1k*Yk/sqrt(scalem), eye(n)] >=  0) ]; 

            %Complementary Sensitivity Constraint
            W2k = W2_w(:,:,j);
            Constraints = [Constraints, ([(ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp)/scalem, ctranspose(W2k*Gk*Xk)/sqrt(scalem); 
                            W2k*Gk*Xk/sqrt(scalem), eye(n)] >= 0 )];        
            %Input Sensitivity Constraint
            W3k = W3_w(:,:,j);
            Constraints = [Constraints, ([(ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp)/scalem, ctranspose(W3k*Xk)/sqrt(scalem); 
                            W3k*Xk/sqrt(scalem), eye(n)] >= 0 )];   
        end        
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLVING THE PROBLEM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % Set some options for YALMIP and solver
    options = sdpsettings('verbose',2,'solver','mosek'); 
    options.mosek.MSK_DPAR_DATA_TOL_AIJ = 1e-16;
    % Solve the problem
    sol = optimize(Constraints,Objective,options)
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANALYZE SOLUTION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
    if (sol.problem == 1 || sol.problem == 9)
     display('Problem Infeasible');
     sol.info
     yalmiperror(sol.problem)
     break;
    end
    
    Yinv = tf(zeros(n,n));
    for j=1:(ny+1)
        Yinv = Yinv + value(Y(:,:,j))*z^(j-1);
    end
    Xstab = tf(zeros(n,n));
    for j=1:(nx+1)
        Xstab = Xstab + value(X(:,:,j))*z^(j-1);
    end
    if(max(abs(zero(Yinv))) > 1 || max(abs(zero(Xstab))) > 1 )
       display('Warning - Unstable Pole or Zero appeared in Controller!')
        break;
    end        
            
    display(k) %to see progress...
    if abs(obj-value(Objective))/value(Objective)<0.05 % set convergence criterion
        display('Converged!')
        obj = value(Objective)  
        break;
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

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL CONTROLLER & ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

X = tf(zeros(n,n)); Yinv = tf(zeros(n,n));
for j=1:(ny+1)
    Yinv = Yinv + Y_c(:,:,j)*z^(j-1);
end
Yinv = inv(Yinv);
% Yinv = inv(Yinv.*Fy);

for j=1:(nx+1)
    X = X + X_c(:,:,j)*z^(j-1);
end
X = X.*Fx;
K = X*Yinv;

status = '.';
if sol.problem == 0 || 4
    status = 'Converged!';
elseif sol.problem == 1
    status = 'Infeasible!';
elseif sol.problem == 9
    status = 'Solver Error! Is the License File of Mosek expired?';
end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE CONSTRAINTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% G_d = frd(G{1}); 
% G_d.Ts = Ts;
% Kfrd = K/(z-1);
% 
% % Characteristic TFs
% S={};T={};U={};
% for i=1:nmod
%     S{i}=feedback(eye(n),G_d*Kfrd);
%     T{i}=feedback(G_d*Kfrd,eye(n));
%     U{i}=feedback(Kfrd,G_d);
% end
% 
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
%     sigma(inv(W3),'y--')
%     title('U')
%     
%     figure(4)
%     bodemag(G_d*Kfrd);
%     hold on
%     bodemag(Ld,wobj)
% end

end