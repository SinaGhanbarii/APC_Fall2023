% Assume All the reactions are 1st order
close, clear, clc

kf = [0.1 0.01 0.05 0.06 0.02];
Vol = 10;
n0 = [110, 1, 0, 2];
C0 = n0/Vol;
tspan = 0:0.1:1000;
%CA = ones(1,length(t));
%CB = ones(1,length(t));
%CC = ones(1,length(t));
%CD = ones(1,length(t));

%ODE Resolution
[t,C] = ode45(@(t,C)ODEsystem(t,C,kf),tspan, C0);

% Primary Results
CA = C(:,1);
CB = C(:,2);
CC = C(:,3);
CD = C(:,4);

% Analyzing Results
CBmax = max(CB);
indx = find(CB == CBmax);
t_max = t(indx);
Conversion_A = (C0(1)-CA)./C0(1) * 100;
selectivity_B = (CB-C0(2))./(C0(1)-CA)*100;
yield_B = CB/C0(1) * 100;
Conv_CBmax = Conversion_A(indx);
Sel_CBCmax = selectivity_B(indx);
Yield_CBmax = yield_B(indx);
disp('the maximum yield of B is:')
disp(Yield_CBmax)

% Visualization
figure
plot(t,C)
legend('C_A','C_B','C_C','C_D')
xlabel('time (s)')
ylabel('Concentration (mol/L)')

figure(2)
plot(t,Conversion_A)
xlabel('time(s)')
ylabel('Conversion of A (%)')


%% Function
function dC = ODEsystem(t,C,k)
    CA = C(1);
    CB = C(2);
    CC = C(3);
    CD = C(4);

    r1 = k(1)*CA;
    r2 = k(2)*CB;
    r3 = k(3)*CB;
    r4 = k(4)*CA;
    r5 = k(5)*CD;

    dCA = -r1+r2-r4+r5;
    dCB = +r1-r2-r3;
    dCC = +r3;
    dCD = r4-r5;

    dC = [dCA dCB dCC dCD]';

end