function [sys, xinit, str, ts] = heatedTankFcn(t, x, u, flag)

% task1 - nominal Qel = 60 W
% Steady State calculation - thetaL_s, thetaH_s 

% task2 - required thetaL = 65°C
% Steady State calculation - Qel_s, thetaH_s 
% Recalculate initial conditions with Qel_s - thetaL_s, thetaH_s 

% task3 - step change (from task1 steady state "initial conditions" 
%                        to task2 steady state "new steady state")
%% template S-function model 

% inputs  - Qel,thetaE
% states  - thetaL(t),thetaH(t)
% outputs - thetaL(t),thetaH(t)

%% parameters - system dimensions, model parameters

nu = 2; % number of inputs
nx = 2; % number of states
ny = 2; % number of outputs

% model parameters
Qel_s = 60; %W
thetaE_s = 20; %degrees celsius

rho_L = 1000; %kg/m3
cp_L = 4180; %J/kg/K
V_L = 0.01; %m3
m_L = rho_L*V_L; %kg
m_H = 1; %kg
A_H = 0.03; %m2
cp_H = 466; %J/kg/K
alpha_LH = 1300; %W/m2/K
alpha_LE = 25; %W/m2/K
a = V_L^(1/3); %m
A_L = a^2; %m2

%time constants
T_H = (m_H*cp_H)/(alpha_LH*A_H);
T_L = (m_L*cp_L)/(alpha_LH*A_H + alpha_LE*A_L);

%gains
K_elH = 1/(alpha_LH*A_H);
K_HL = (alpha_LH*A_H)/(alpha_LH*A_H + alpha_LE*A_L);
K_EL = (alpha_LE*A_L)/(alpha_LH*A_H + alpha_LE*A_L);

thetaL_s = 65; %degrees celsius
Qel2_s = 52.2179; %W
thetaL_ref  = 75;

Kc1         = 1;
Kc2         = 10;
Kc3         = 100;
%% initial conditions x(0) = x0
A = [-1 1;
     1 -K_HL];
b = [K_elH*Qel_s;
     K_EL*thetaE_s];

%task2
%A2 = [-K_elH 1;
 %     0 -K_HL];
%b2 = [thetaL_s ;
%      K_EL*thetaE_s - thetaL_s];

% task 4
A4 = [1 -K_HL;
    Kc2*K_elH-1 1];
b4 = [K_EL*thetaE_s;
    K_elH*(Qel_s+Kc2*thetaL_ref)];

x0  = A\b; % vector of initial conditions
%x02 = A2\b2;
%x02 = A4\b4;
%% dynamics equations dx(t)/dt = f(x(t), u(t))

    function f = dynamics_equations (x, u)
        thetaL = x(1); 
        thetaH = x(2);
        Qel    = u(1); 
        thetaE = u(2);
        %Task 2:
       % thetaH = x02(1);
       % Qel2 = x02(2);
       % thetaE = u02(1);
       % thetaL = u02(2);
       
        f = [(K_HL*theta_H + K_EL*theta_E - theta_L)/T_L ;
             (theta_L - theta_H + K_elH*Qel)/T_H];
    end

%% output equations y(t) = g(x(t))

    function g = output_equations (x)
        theta_L = x(1);
        theta_H = x(2);
        
        g = [theta_L ; theta_H]';
    end

%% END OF S-FUNCTION














%% other functions // do not change unless really needed
switch flag
    case 0
        [sys, xinit, str, ts] = mdlInitializeSizes;
    case 1
        sys = mdlDerivatives(t,x,u);
    case 3
        sys = mdlOutputs(t,x,u);
    case {2, 4, 9}
        sys = [];
end

    function [sys, xinit, str, ts] = mdlInitializeSizes
        sizes = simsizes;
        sizes.NumSampleTimes = 1;
        sizes.NumContStates  = nx; % pocet spojitych stavov
        sizes.NumOutputs     = ny; % pocet vystupov
        sizes.NumInputs      = nu; % pocet vstupov
        xinit = x0; % zaciatocne podmienky pre dif. rovnice
        sys = simsizes(sizes);
        str = []; ts = [0 0];
    end
    function f = mdlDerivatives(~, x, u)
        f = dynamics_equations (x, u);
    end
    function y = mdlOutputs(~, x, ~)
        y = output_equations(x);
    end
end
%%
