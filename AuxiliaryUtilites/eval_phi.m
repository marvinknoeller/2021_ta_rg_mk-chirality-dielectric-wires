function [chir,smooth_relax,val,FFM,dist,pen1,pen2] = eval_phi(Var,x,t,d)
    n = length(x);
    if nargin == 2
        curve = x;
    elseif nargin ==4
        curve = x+t*d;
    else
        error("function is only defined for 2 or 4 inputs")
    end
    FFM = FarFieldMatrixFunction_Spline(curve,Var);
    [chir,smooth_relax] = chiral(FFM);
    Len = sum(sqrt(sum( (curve(:,2:end)-curve(:,1:end-1)).^2,1 )));
    dist = Len/(n-1);
    pen1 = Seg_Pen2(curve,Var.M,n,Var.length);
    pen2 = Curvature_Constraint(curve,Var.M,n);
    val = -(smooth_relax - Var.lambda1*pen1 - Var.lambda2*pen2);
end