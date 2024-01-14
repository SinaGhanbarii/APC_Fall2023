clear; close all; clc;

% data
p = 1.013;                  % [bar]
R = 8.314*10^-2;

x = [0.047  0.107  0.845];  % [MeOH Benz MeCyP]

r  = [2.11  3.19  3.97];
q  = [1.19  2.40  3.01];

a_ij = [-128.9  -118.3  -6.47 ];    % [12 - 13 - 23]
a_ji = [ 997.4   1384    56.47];    % [21 - 31 - 32]

% create a matrix to condensate the formulas
a = [0 a_ij(1) a_ij(2) ;
    a_ij(1) 0 a_ji(2);
    a_ji(2) a_ji(3) 0];


% fsolve
FG = [340 0.2 0.1 0.6];
sol = fsolve(@(u)uniq(u,a,r,q,x,p,R),FG);
T = sol(1);
y = sol(2:4);

disp(T)
disp(y)
% function
function f = uniq(u,a,r,q,x,p,R)
T = u(1);
y(1) = u(2);
y(2) = u(3);
y(3) = u(4);
% data
B = [-1155 -587  -618;
     -587  -1086 -1134;
     -618  -1134 -1186]./10^3;

v_L = [61.1  93.7  118]./10^3;

% virial EoS
B_mix = B(1,1)*y(1)^2+B(2,2)*y(2)^2+ B(3,3)*y(3)^2+ 2*B(1,2)*y(1)*y(2) + ...
    2*B(1,3)*y(1)*y(3)+ 2*B(2,3)*y(2)*y(3);

% writing with a for cycle
ln_phi = ones(1,3);
for i = 1:3
    ln_phi(i) = p/(R*T) * (2*y(1).*B(i,1)+y(2).*B(i,2)+y(3).*B(i,3) - B_mix);
end

phi_V = exp(ln_phi);

% data antoine
Ant_A = [4.92531   4.01814	 3.98773 ];
Ant_B = [1432.526  1203.835  1186.059];
Ant_C = [-61.819   -53.226	 -47.108 ];

P_ev = 10.^(Ant_A-(Ant_B./(Ant_C+T)));

% pure liquid
phi_LS   = exp(p/(R*T) .* [B(1,1) B(2,2) B(3,3)]);
poynting = exp(v_L.*(p-P_ev./(R*T)));
fug_L    = P_ev.*phi_LS.*poynting;


% uniquac (vapor)
z = 10;

PHI   = x.*r ./ (x(1)*r(1) + x(2)*r(2) + x(3)*r(3));
THETA = x.*q ./ (x(1)*q(1) + x(2)*q(2) + x(3)*q(3));

L = z*(r-q)/2 - (r-1);

TAU = exp(-a./T);

ad_1 = log(PHI./x);
ad_2 = z.*q.*log(THETA./PHI)/2;
ad_3 = L;
ad_4 = -(PHI./x)*sum(x.*L);
ad_6 = q;

% extended writing with for cycle
ad_5 = ones(1,3);
for i = 1:3
    ad_5(i) = -q(i)* log(THETA(1)*TAU(1,i) + THETA(2)*TAU(2,i) + THETA(3)*TAU(3,i));
end

ad_7 = ones(1,3);
for i=1:3
    ad_7(i) = -q(i)*((THETA(1)*TAU(i,1))/ (THETA(1)*TAU(1,1) + THETA(2)*TAU(2,1) + THETA(3)*TAU(3,1) ...
        )) + (THETA(2)*TAU(i,1))/ (THETA(1)*TAU(1,2) + THETA(2)*TAU(2,2) + THETA(3)*TAU(3,2)) ...
        + (THETA(3)*TAU(i,1))/ (THETA(1)*TAU(1,3) + THETA(2)*TAU(2,3) + THETA(3)*TAU(3,3));
end

ln_gamma = ad_1 + ad_2 + ad_3 + ad_4 + ad_5 + ad_6 + ad_7;
gamma    = exp(ln_gamma);

f(1) = x(1)*gamma(1).*fug_L(1) - (y(1).*phi_V(1).*p);
f(2) = x(2)*gamma(2).*fug_L(2) - (y(2).*phi_V(2).*p);
f(3) = x(3)*gamma(3).*fug_L(3) - (y(3).*phi_V(3).*p);

f(4) = y(1) + y(2) + y(3) - 1;
end