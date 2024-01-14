% Practical Session 1 - Excercise 2 - Flash Evaporation
close, clear, clc
options = optimoptions('fsolve','MaxIter',1e10,'MaxFunEval',1e10,'Functiontolerance',10^-10,'display','off');
F  = 10;     % [mol/s]
z  = [0.25  0.25   0.3    0.2];
T0 = 500;   % [K]
P0 = 4;     % [bar]
Pf = 3;
Tf = 400;

% log10(Pi0(T)) = A - B /(T+C)
% species: [-  Hexane  Octane  Decane]  =  [- 2 3 4]
A = [4.00266 4.04867 4.07857];
B = [1171.530 1355.126 1501.268];
C = [-48.784 -63.633 -78.670];

z_uncond = z(1);
z_cond   = z(2:4);

alpha_DP = 1;
alpha_BP = z_uncond;

% fsolve of Rachford-Rice Equation
Td = fsolve(@(T)UncondRR(alpha_DP,z_uncond,z_cond,T,P0,A,B,C),500,options);
Tb = fsolve(@(T)UncondRR(alpha_BP,z_uncond,z_cond,T,P0,A,B,C),500,options);
alpha = fsolve(@(alpha_DP)UncondRR(alpha_DP,z_uncond,z_cond,Tf,Pf,A,B,C),1,options); %alpha = V/F
z_out = fsolve(@(z_uncond)UncondRR(alpha_DP,z_uncond,z_cond,Tf,Pf,A,B,C),z(1),options); %alpha = V/F

V = alpha*F;
L = (1-alpha)*F;
disp('Td is: ')
disp(Td)
disp('')
disp('Tb is: ')
disp(Tb)
disp('')
disp('The Output Liquid and Vapor is: ')
disp([L,V])
result = newRR(alpha,z_uncond,z_cond,Tf,Pf,A,B,C);
disp(result)
function F = UncondRR(alpha,z_uncond,z_cond,T,P0,A,B,C)
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
    F = sum(x-y);
end

