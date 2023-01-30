function [der_f] = eval_Gradf(Var,x,adder,FarFieldMatrix,t,d)
    n = length(x);
    if nargin == 4
        curve = x;
    elseif nargin == 6
        curve = x+t*d;
        [~,~,~,FarFieldMatrix] = eval_phi(Var,x,t,d);
    else
        error("function is only defined for 2 or 4 inputs")
    end
    full_der = zeros(3,length(curve));
    parfor kk = 1:length(curve) 
        for j = 1:3
            [der_ffm] = derive_farfieldmatrix_Spline(curve,squeeze( adder(:,3*(kk-1)+j,:)),Var);
            [full_der(j,kk)] = derive_measure(FarFieldMatrix,der_ffm);
            J1(j,kk) = Der_Seg_Pen2(curve,squeeze( adder(:,3*(kk-1)+j,:)),Var.M,n,Var.length);
            J2(j,kk) = Der_Curvature_Constraint(curve,squeeze( adder(:,3*(kk-1)+j,:)),Var.M,n);
        end
    end
    
    der_f = -(full_der - Var.lambda1*J1 - Var.lambda2*J2);
end