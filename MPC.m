%% ================ Project 2 ====================
clear;clc
%Initialization
load("New_demand.mat");
demandM = [];
for i = 1:length(ts_power.Data)
    demandM(i) = ts_power.Data(i);
end

load("New_generation.mat");

pvM = [];
for i = 1:length(ts_power.Data)
    pvM(i) = ts_power.Data(i);
end

d = pvM - demandM;
A = eye(2);
B = 1e-3*[1.56 1.56; -5.66 0];
Bd = [1.56e-3; 0];
C = eye(2);
gamma = [1e-8 1e-8];         % - Weighting faactor setpoint tracking
lambda = [4 1e-3];    % - Weighting factor of efforts to increase control
alpha = [5e-3 8e-3];       % - Weighting factor control effort PH2, Pgrid
NT = length(pvM); Np = 10; Nc = 2;
nx = size(A,1); % number of states
nu = size(B,2);% number of control inputs
ny = size(C,1);% number of outputs
nd = size(Bd,2);

%% Convert to incremental model
A_in = [A B; zeros(nu,nx) eye(nu,nu)]; % A_in is M in slides
B_in = [B; eye(nu,nu)]; % B_in is N in slides
Bd_in = [Bd; zeros(nu,nd)];
C_in = [C zeros(nu,nu)]; % C_in is Q in slides
nx_in = size(A_in,1); % number of states of incremental
nu_in = size(B_in,2); % number of inputs
ny_in = size(C_in,1); % number of outputs
nd_in = size(Bd_in,2);

%% Weight matrix of the cost function

gamma1 = [];
for i = 1:1:Np
    gamma1 = [gamma1 gamma];
end
gammaM = diag(gamma1);

lambda1 = [];
for i = 1:1:Nc
    lambda1 = [lambda1 lambda];
end
lambdaM = diag(lambda1);

alpha1 = [];
for i = 1:1:Nc
    alpha1 = [alpha1 alpha];
end
alphaM = diag(alpha1);

%% Prediction
% - Matrix F
F = [];
for i = 1:1:Np
    F = [F; C_in*A_in^(i)];
end
F;

H = [];
zerosH = zeros(size(C_in*B_in));
for i = 1:1:Np
    hLine = [];
    for j = 1:1:Nc
        if j <= i
            eme = A_in^(i-j);
            hLine = [hLine C_in*eme*B_in];
        else
            hLine = [hLine zerosH];
        end
    end
    H = [H;hLine];
end

Hd = [];
zerosHd = zeros(size(C_in*Bd_in));
for i = 1:1:Np
    hdLine = [];
    for j = 1:1:Np
        if j <= i
            eme = A_in^(i-j);
            hdLine = [hdLine C_in*eme*Bd_in];
        else
            hdLine = [hdLine zerosHd];
        end
    end
    Hd = [Hd;hdLine];
end

%% Constraints
Ymax = [90; 70]; Ymin = [10; 40];
YMaxArray = [];
YMinArray = [];
for i = 1:1:Np
    YMaxArray = [YMaxArray;Ymax];
    YMinArray = [YMinArray;Ymin];
end

deltaUmax = [0.05*900; 1*900]; deltaUmin = [-0.05*900; -1*900];
deltaUMaxArray = [];
deltaUMinArray = [];
for i = 1:1:Nc
    deltaUMaxArray = [deltaUMaxArray;deltaUmax];
    deltaUMinArray = [deltaUMinArray;deltaUmin];
end
B1 = deltaUMaxArray;
B2 = -deltaUMinArray;

Umax = 300*[1.2; 6]; Umin = 700*[-2.2; -2.5];
UMaxArray = [];
UMinArray = [];
for i = 1:1:Nc
    UMaxArray = [UMaxArray;Umax];
    UMinArray = [UMinArray;Umin];
end

T=[];
miniT=eye(nu);
miniZeros=zeros(nu);
for i=1:1:Nc
    tLine = [];
    for j=1:1:Nc
        if j<=i
            tLine = [tLine, miniT];
        else
            tLine = [tLine, miniZeros];
        end
    end
    T=[T;tLine];
end

A1 = eye(Nc*nu,Nc*nu);
A2 = -A1;
A3 = H;
A4 = -A3;
A5 = T;
A6 = -A5;

Aqp = [A1;A2;A3;A4;A5;A6]; % Matrix R
%% set point matrix

% Ref matrix
yref = [50; 50];
w = zeros(ny*Np,1);
for i = 1:ny:Np*ny
    w(i:i+ny-1,1) = yref;
end

%% Construction of optimization problem
x0 = [60; 30;0; 0]; % SOC LOH PH2 Pbat
x = zeros(nx_in,NT+1); x(:,1) = x0;
yout = zeros(ny, NT+1); yout(:,1) = C_in*x0;

%%
ut = [100; -200 ]; % PH2, Pbattery

for k = 1: NT -(Np -1)
    utM = ones(nu*Nc,1);
    for j = 1:nu:nu*Nc
        utM(j:j+nu-1,1) = ut;
    end
    xk = x(:,k);
    D = d(k:k+Np-1); D = D';
    
   % J = (y-w)'*gammaM*(y-m) + DeltaU'*lambdaM*DeltaU + U'*AlphaM*U
    Hqp = (H'*gammaM*H + lambdaM + T'*alphaM*T);
    Bqp = (F*xk + Hd*D - w)'*gammaM*H + utM'*alphaM*T;
    
    utArray = [];
    for j = 1:1:Nc
        utArray = [utArray;ut];
    end
    
    B3 = YMaxArray - (F*xk + Hd*D);
    B4 = F*xk + Hd*D - YMinArray;

    B5 = UMaxArray - utArray;
    B6 = utArray - UMinArray;
    
    Brqp = [B1;B2;B3;B4;B5;B6];
    
    [deltaU, fval, exitflag] = quadprog(Hqp,Bqp,Aqp,Brqp,[],[],[],[],[]);
    % u = deltaU(1:nu) + u;
    x(:,k+1) = A_in*xk + B_in*deltaU(1:nu) + Bd_in*D(1);
    yout(:,k+1) = C_in*x(:,k+1);
    
    ut = ut +  deltaU(1:nu,1); 
end
%%
time = (0:NT-1); figure(1); hold on; 
plot(pvM,'r','LineWidth',.7);
plot(demandM, 'b','LineWidth',.7);
plot(time, x(3,1:86400),'g','LineWidth',.7)
plot(time, x(4,1:86400),'k','LineWidth',.7)
legend('Generation','Demand', 'PH2', 'Pbattery')
xlabel('Time(s)'); ylabel('Power(W)');
title("Microgrid System Power Flow");
%%
time2 = (0:NT);
figure(2); hold on;
plot(time2, yout(1,1:86401), 'r', 'LineWidth', .7)
plot(time2, yout(2,1:86401), 'b', 'LineWidth', .7)
legend('LOH','SOC')
xlabel('Time(s)'); ylabel("%");
title("Storage devices");
