function name_of_file = BFGS_Backtrack(example,N,h,curve_length,save_results)
folder = fileparts(which(mfilename));
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
%% technical setup -- just for programming
% N = 5;
for nn=1:N
    NVec(nn^2:nn^2+2*nn)=nn;
end
% size of vectors and matrix
Q = 2*N*(N+2);
Qd2 = Q/2;
%% Initial guess
subnum = 19;    %number of subsegments
num_x = subnum+1;   %number of points that define my curve
t = linspace(0,1,num_x);

% example = 1;
switch example
    case 1
        %% EXAMPLE 1
        rng(123) %for different orientation use 123456
        eps_rel = 5;
        mu_rel = 1;
        kappa = 1;
        roh = 0.03;
        %         h = 6;
        R = 0.;
        x0 = R*cos(0.5*pi*t)+ 0.02*(2*rand(1,num_x)-1);
        y0 = R*sin(0.5*pi*t)+ 0.02*(2*rand(1,num_x)-1);
        z0 = h*t-h/2;
        p = [x0;y0;z0];
        plot3(p(1,:),p(2,:),p(3,:),'-b','LineWidth',5);
        lambda1 = 0.5;
        lambda2 = 5e-4;
        Var.length = curve_length;
    case 2
        %% EXAMPLE 2
        rng(123)
        eps_rel = 30;
        mu_rel = 1;
        kappa = 1;
        roh = 0.03;
        %         h = 10;
        R = 1;
        x0 = -R * ones(1,num_x) + 0.02*(2*rand(1,num_x)-1);
        y0 = 0 * zeros(1,num_x) + 0.02*(2*rand(1,num_x)-1);
        z0 = -(h*t -h/2);
        p1 = [x0;y0;z0];
        
        x0 = R * ones(1,num_x) + 0.02*(2*rand(1,num_x)-1);
        y0 = 0 * zeros(1,num_x) + 0.02*(2*rand(1,num_x)-1);
        z0 = (h*t -h/2);
        
        p2 = [x0;y0;z0];
        t = linspace(0.2,0.8,5);
        x0 = R*cos(pi*t+pi)+ 0.02*(2*rand(1,5)-1);
        y0 = 0*zeros(1,5)+ 0.02*(2*rand(1,5)-1);
        z0 = R*sin(pi*t+pi)-h/2;
        
        bow = [x0;y0;z0];
        p = [p1,bow,p2];
        lambda1 = 0.1;
        lambda2 = 1e-4;
        Var.length = curve_length;
end

%% Var
Var.kappa = kappa;
Var.M = 5;
Var.eps_rel = eps_rel;
Var.mu_rel = mu_rel;
Var.roh = roh;
Var.Q = Q;
Var.NVec = NVec;
n = length(p);
adder = zeros(3,3*n,n);
[P, tspline, coefs, br,ts] = splinepoints(p,Var.M);
plot3(P(1,:), P(2,:) ,P(3,:),'-r')
for kk = 1:n
    adder(1:3,3*(kk-1)+1:3*kk,kk) = 1*eye(3);
end
B_Hess = eye(3*n);

pp(:,:,1) = p;

grad_f(:,:,1) = zeros(3,n);
Var.lambda1 = lambda1;
Var.lambda2 = lambda2;
num_x = size(p,2);
%%
string_length = num2str(Var.length);
string_point = string_length == '.';
string_length(string_point) = '_';
name_of_file = strcat(num2str(example),"_EX",string_length,"_L");
%%
for steps = 1 : 500
    %% Evaluate to get the FFM
    [chir(steps),smooth_relax(steps),func,FarFieldMatrix,~,pen1,pen2] = eval_phi(Var,p);
    disp(strcat('Chiralitaetsmass: ' ,{' '}, num2str(chir(steps))));
    disp(strcat( 'Smooth Relaxation: ',{' '}, num2str((smooth_relax(steps)))));
    [~,~,cint] = chiral(FarFieldMatrix);
    [~,Len] = CurvePenalty(p,Var.M,num_x);
    disp(strcat( 'length of curve: ',{' '}, num2str((Len))));
    %% Store interesting values
    ValVec(steps) = smooth_relax(steps);
    Pen1Vec(steps) = lambda1*pen1;
    Pen2Vec(steps) = lambda2*pen2;
    CintVec(steps) = cint;
    %%
    if steps == 1
        Gradf = eval_Gradf(Var,p,adder,FarFieldMatrix);
        grad_f(:,:,steps) = Gradf;
    else
        Gradf = grad_f(:,:,steps);
    end
    % Step 1 : Solve the linear system
    jac_shaped = reshape(Gradf,[],1);
    p_k_shaped = linsolve(B_Hess,-jac_shaped);
    p_k = reshape(p_k_shaped,3,n);
    
    % Initialization
    alpha = 1e-4;
    phi_0 = func;
    Dphi_0 = jac_shaped.'*p_k_shaped;
    jj = 0;
    lambda = 0.9;
    [~,~,phi_val] = eval_phi(Var,p,lambda^jj,p_k);
    while phi_val > func + alpha * lambda^jj * Dphi_0
        jj = jj+1;
        [~,~,phi_val] = eval_phi(Var,p,lambda^jj,p_k);
    end
    lambda_fin = lambda^jj;
    X1 = p + lambda_fin*p_k;
    total_movement = sum(sqrt((sum(p - X1).^2)));
    display([steps, total_movement])
    p = X1;
    FarFieldMatrix = FarFieldMatrixFunction_Spline(p,Var);
    grad_f(:,:,steps+1) = eval_Gradf(Var,p,adder,FarFieldMatrix);
    pp(:,:,steps + 1) = p;
    s_k = (pp(:,:,steps+1) - pp(:,:,steps));
    s_k_shaped = reshape(s_k,[],1);
    s_k_shapedT = s_k_shaped.';
    if norm(s_k_shaped)/norm(X1-lambda_fin*p_k)<1e-4
        fprintf("too little movement. Break.\n");
        break
    end
    % 4th Step: y_k
    gamma_k = (grad_f(:,:,steps+1) - grad_f(:,:,steps));
    gamma_k_shaped = reshape(gamma_k,[],1);
    gamma_k_shapedT = gamma_k_shaped.';
    norm_g = sqrt(sum(sum(Gradf.^2,1)));
    
    y_k = (grad_f(:,:,steps+1) - grad_f(:,:,steps));
    y_k_shaped = reshape(y_k,[],1);
    y_k_shapedT = y_k_shaped.';
    if steps == 1
        init_scal = (y_k_shapedT*s_k_shaped)/(s_k_shapedT*s_k_shaped);
        fprintf(strcat("initial scaling size of the approximate Hessian: ", num2str(init_scal),"\n"));
    end
    if s_k_shapedT*y_k_shaped<0
        fprintf("Might destroy spd\n")
    end
    % 5th Step: update the matrix
    if (s_k_shapedT*y_k_shaped/(s_k_shapedT*s_k_shaped)) >= 1e-5*norm_g
        B_Hess = B_Hess + y_k_shaped*y_k_shapedT/(y_k_shapedT*s_k_shaped) ...
            - B_Hess*s_k_shaped*s_k_shapedT*B_Hess/(s_k_shapedT*B_Hess*s_k_shaped);
    else
        fprintf("No Update\n");
        B_Hess = B_Hess;
    end
    if min(eig(B_Hess))<0
        fprintf("Matrix not spd\n")
    end
    p_spline = splinepoints(p,Var.M);
    subplot(1,2,1)
    plot3(p_spline(1,:),p_spline(2,:),p_spline(3,:),'-*','LineWidth',2)
    subplot(1,2,2)
    plot3(p_spline(1,:),p_spline(2,:),p_spline(3,:),'-*','LineWidth',2)
    view(2)
    axis square
    drawnow
end

figure
plot(ValVec,'rd')
hold on
plot(Pen1Vec,'go')
plot(Pen2Vec,'b*')

if save_results == 1
    save(name_of_file);
end
end
