clear; close all; clc;
format long

% miscela acetone (4.7 mol%) e n-pentano (95.3%)
% Acetone:  CH3 - CH3CO     [a - b] = [1(1) 9(18)] [main(sub)]
% n-Petano: CH3 - CH2       [c - d] = [1(1) 1(2)]  [main(sub)]

% data
T  = 307;       % [K]
p  = 1.013;     % [bar]

x1 = 0.047;     
x2 = 1-x1;
x  = [x1 x2];

z = 10;

R1  = 0.9011;  
R18 = 1.6724;
R2  = 0.6744;

Q1  = 0.848;
Q18 = 1.488;
Q2  = 0.540;

% number of groups in each molecule
nu_1 = ;
nu_2 = ;

% order of the groups in the molecules
R_1  = ;   
R_2  = ;

Q_1  = ;   
Q_2  = ;


% preliminary calculations (here we just consider the molecules, not the groups!!)
r1 = ;      
r2 = ;
r  = [r1 r2];

q1 = ;      
q2 = ;
q  = [q1 q2];

Phi1 = ;   
Phi2 = ;
Phi  = [Phi1 Phi2];

Theta1 = ; 
Theta2 = ; 
Theta  = [Theta1 Theta2];

L1 = ;  
L2 = ;
L  = [L1 L2];

% combinatorial contribution to gamma (here we have just considered the MOLECULES!)
ln_gammaC = log(Phi./x) + 0.5*z.*q.*log(Theta./Phi) + L - Phi./x .*sum(x.*L);




% % residual contribution to gamma (here we start considering the GROUPS!!)
% data
a_19   = 476.40;    % a_ch3-ch3co           
a_91   = 26.760;    % a_ch3co-ch3


% in psi evaluation, we just consider the differences in the main groups!
psi_19 = ;      
psi_91 = ;

psi = [1        1       psi_91;
       1        1       psi_91;
       psi_19   psi_19  1];


% residual activity coefficient in solution with only 1 molecular species
% solution with just acetone (1)
X1_1   = ;
X1_18  = ;
X1     = [X1_1  X1_18];

THETA1_1  = ;
THETA1_18 = ;
THETA1    = [THETA1_1  THETA1_18];

% extended writing
lnDELTA1_1 = ;

lnDELTA1_18 = ;

lnDELTA1    = [lnDELTA1_1  lnDELTA1_18];


% solution with just pentane (2)
X2_1 = ;
X2_2 = ;
X2   = [X2_1  X2_2];
% since in n-pentane there's only the MAIN group n.1 (CH3) we have:
lnDELTA2_1 = ;
lnDELTA2_2 = ;
lnDELTA2   = [lnDELTA2_1  lnDELTA2_2];


% back to the mixture (acetone + pentane)
% compute the molar fractionof the groups in the mixture:
X_g1  = ;
X_g2  = ;
X_g18 = ;

theta_g1  = ;
theta_g2  = ;
theta_g18 = ;

% order the values of theta_gi and psi_ij in order to write the summations more easily
theta_g = [theta_g1 theta_g2 theta_g18];

psi_r1  = psi(1,:);
psi_r2  = psi(2,:);
psi_r3  = psi(3,:);

% condensed writing
lnDelta_g1  = ;

lnDelta_g2  = ;

lnDelta_g18 = ;


% final calculation of the residual activity coefficients
ln_gamma1_R = ;
ln_gamma2_R = ;
ln_gammaR   = [ln_gamma1_R   ln_gamma2_R];


% final calculation of the total activity coefficients
ln_gamma = ;
gamma = ;




