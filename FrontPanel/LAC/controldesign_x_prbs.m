function [K,Ld,status] = controldesign_x_prbs(G_freqresp, Ts, fc)

s = tf('s');
z = tf('z',Ts);

%% Plant
wobj = logspace(log10(100*2*pi),log10(1/Ts*pi),500);

wmin = wobj(1);
wmax = wobj(end);

G{1} = G_freqresp / abs(freqresp(G_freqresp,wobj(1)));
G{1}.Ts=Ts;
nmod=1; %number of models
n = size(G{1},1);

%% initial controller 
if real(freqresp(G{1},wobj(1))) < 0
    Kinit = -tf(1e-3);
else
    Kinit = tf(1e-3);
end

%% Controller Structure and Initial Controller
nx = 10; ny = 10; %Controller Order

for j=1:n
    Kden = [pole(Kinit(j,j))', 0*ones(1, ny+1-length(Kinit(j,j).den{1}))] ;
    Y_c(j,j,:) = flip(poly(Kden));
    Knum = flip(Kinit(j,j).num{1} * poly(0*ones(1,nx))) ;
    X_c(j,j,:) = dcgain(Kinit(j,j))/(sum(Knum)/sum(poly(Kden)))*Knum;   
end

% Fx, Fy are multiplied element-wise!
Fy = tf(ones(n));%*(z-1)*tf(a,1,Ts);
Fx = tf(ones(n));%*tf(b,1,Ts);

% load K_cl
% X_c = X_cl{end}; Y_c = Y_cl{end};

%% Control Performance

%Desired open-loop transfer function Ld
[b,a] = butter(2, fc*2*Ts);
Ld = tf(b,a,Ts);

% Open-loop TF constraint
clear Wol
tmp = freqresp(Ld,wobj);
tmp(end-4:end) = tmp(end-5);
for i=1:length(wobj)
    Wol(:,:,i) = inv( tmp(i)*diag(ones(1,n)) );
end
Wol = frd(Wol,wobj);

%% Generate Frequency Responses
for i=1:nmod
    G_w{i} = freqresp(G{i},wobj);
end
G_w{end}(end)=G_w{end}(end-1);

Ld_w = freqresp(Ld,wobj);
Wol_w = freqresp(Wol,wobj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Iterative Algorithm
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
    % Minimize Hinf norm between L and Ld
    Objective=0;
    %Hinf
    gamma_LS = sdpvar(n,n, nmod*length(wobj),'hermitian');
    for j=1:nmod*length(wobj)
        Objective = Objective + trace(gamma_LS(:,:,j)); %maximize gamma => minimize -gamma   wobj(j)/wobj(end)*
    end      
    Constraints = [Constraints, gamma_LS >= 0];        
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRIDDING CONSTRAINTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for jj=1:nmod       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Performance Constraints    
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

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % Stability Constraint
            tmp2 = ctranspose(Yk_c_stab)*Yk_c_stab; % Make sure tmp = Yk_c^*Yk_c really is hermitian
            tmp = diag(real(diag(tmp2)));
            tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));
            
            Constraints = [Constraints, [ctranspose(Yk_stab)*Yk_c_stab+ctranspose(Yk_c_stab)*Yk_stab-tmp] >= 1e-10 ];            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            %Objective Function
            scalem=1;
            %Minimize ||L-Ld||_2
            Ldk = Ld_w(:,:,j);  
            tmp2 = ctranspose(Yk_c)*Yk_c;
            tmp = diag(real(diag(tmp2)));
            tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));
            scalem = 1;%abs(tmp);
            Constraints = [Constraints, [ctranspose(Yk)*Yk_c+ctranspose(Yk_c)*Yk-tmp, ctranspose(Gk*Xk-Ldk*Yk); ...
                            Gk*Xk-Ldk*Yk, gamma_LS(:,:,j + (jj-1)*length(wobj))] >=0 ];
                        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            %Open-loop TF constraint
            if wobj(j) >= fc*2*pi
                Wolk = Wol_w(:,:,j);
                Constraints = [Constraints, ([(ctranspose(Yk)*Yk_c + ctranspose(Yk_c)*Yk - tmp)/scalem, ctranspose(Wolk*Gk*Xk)/sqrt(scalem); ...
                                Wolk*Gk*Xk/sqrt(scalem), eye(n)] >=  0)]; 
            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            % Additional stability constraint
            Pk = Yk + Gk*Xk;
            Pk_c = Yk_c + Gk*Xk_c;
            tmp2 = ctranspose(Pk_c)*Pk_c; %Make sure tmp = Pk_c^*Pk_c really is hermitian
            tmp = diag(real(diag(tmp2)));
            tmp = tmp + triu(tmp2,1) + ctranspose(triu(tmp2,1));
          
            Constraints = [Constraints, [ctranspose(Pk)*Pk_c + ctranspose(Pk_c)*Pk - tmp] >= 1e-10 ];                        
        end        
    end
    
     % Set some options for YALMIP and solver
    options = sdpsettings('verbose',2,'solver','mosek'); % mosek sedumi
    % Solve the problem
    sol = optimize(Constraints,Objective,options)
    
    % Analyze error flags
    if sol.problem == 1 || sol.problem == 9
%      sol.info
%      yalmiperror(sol.problem)
     break;
    end
    
    Yinv = tf(zeros(n,n));
    for j=1:(ny+1)
        Yinv = Yinv + value(Y(:,:,j))*z^(j-1);
    end
    if(max(abs(zero(Yinv))) > 1 )
    %    display('Warning - Unstable Pole appeared in Controller!')
        break;
    end      
    
    display(k) %to see progress...
    if abs(obj-value(Objective))/value(Objective)<0.025 % set convergence criterion
        value(Objective)
        display('Converged!')
        break; %convergence!
    else
        obj = value(Objective)        
        X_c = value(X);
        Y_c = value(Y);
        X_cl{k} = value(X);
        Y_cl{k} = value(Y);
        save K_cl X_cl Y_cl       
    end
end
toc

%% Final Controller & Analysis
X = tf(zeros(n,n)); Yinv = tf(zeros(n,n));
for j=1:(ny+1)
    Yinv = Yinv + Y_c(:,:,j)*z^(j-1);
end
Yinv = inv(Yinv);
% Yinv = inv(Yinv.*Fy);

for j=1:(nx+1)
    X = X + X_c(:,:,j)*z^(j-1);
end
% X = X.*Fx;
K = X*Yinv;%.*[1/Fy 1; 1 1/Fy]; 
K = K/dcgain(K);

status = '.';
if sol.problem == 0 || 4
    status = 'Converged!';
elseif sol.problem == 1
    status = 'Infeasible!';
elseif sol.problem == 9
    status = 'Solver Error! Is the License File of Mosek expired?';
end


end