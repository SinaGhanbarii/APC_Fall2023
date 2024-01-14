function m = newRR(alpha,z_uncond,z_cond,T,P0,A,B,C)
    Pi0_T = 10.^(A-B./(T+C));
    ki = Pi0_T/P0;

    %Liquid
    x_cond = z_cond./(1+alpha.*(ki-1));
    x_uncond = 0;
    x = [x_uncond x_cond];

    %Vapor
    y_cond = x_cond.*ki;
    y_uncond = z_uncond./alpha;
    y = [y_uncond y_cond];

    m = [x;y];
end
