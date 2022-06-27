%% Parameter settings for simulation
% Please see the method section: "Computational modeling of ROS trigger
% waves" for parameter details.

% Diffusion coefficient (D): D
% Intercellular distance: h
% Positive feedback strength (kpositivefb): k_pfb
% EC50 of positive feedback (EC50positivefb): EC50_pfb
% Hill coefficient of positive feedback (npositivefb): n_pfb
% Double-negative feedback strength (kdnegativefb): k_dnfb
% EC50 of Double-negative feedback (EC50dnegativefb): EC50_dnfb
% Hill coefficient of Double-negative feedback (ndnegativefb): n_dnfb
% basal production of GSH (cGSH): c_GSH
% EC50 of erastin (EC50erastin): EC50_erastin
% Erasin concentration (E): E
% ROS degradation rate (kdeg): k_deg
% ROS basal production rate (ksynth): k_synth
% Cell number per each side of a square domain: nCell
% Initial and final time of simulation: tStart and tFinal
% Time delay for setting kpositvefb, kdnegativefb, and ksynth to be zero: delay
parameters = struct('D',178,'h',5,'k_pfb',1.2,'EC50_pfb',1,'n_pfb',3,...
    'k_dnfb',1.5,'EC50_dnfb',2,'n_dnfb',3,...
    'c_GSH',0.1,'EC50_erastin',0.27,'E',10,...
    'k_deg',0.26,'k_synth',0.1,...
    'nCell',200,'tStart',0,'tFinal',300,'delay',30);

% Calculate steady states (SS) of ROS
SS = steadystates(parameters); 

% Set initial values of all cells (ros0) to be the lower ROS steady state
parameters.ros0 = SS(1)*ones(parameters.nCell^2,1);

% Set the ROS level within the photoinduction area to be a constant ROS
% (2.57), which is the mean of USS and HSS when erastin = 10 µM
mask = zeros(parameters.nCell);
mask(1:12,:) = 1; % for the radius of 248 µm, 12 cells on the radius
parameters.idx_center = find(mask);
parameters.ros0(parameters.idx_center) = 2.57;

% Set the death threshold to be 90% of the higher steady state if there are
% 3 steady-states; otherwise set the death threshold to be 4.26, which is
% 90% of the higher steady state when erastin = 10 µM.
if length(SS)==1 % for monostable conditon
    parameters.death_threshold = 4.26; 
else % for bistable conditon
    parameters.death_threshold = 0.9*SS(3);
end

% Set boundary condition of simulation to be zero
parameters.ROS_boundary = 0*[1 1 1 1];

%% Wave simulation
fcn_1([],[]) % record cell death events
options = odeset('Refine',4,'NonNegative',1:length(parameters.ros0));
sol = ode45(@(t,y) ftw2D(t,y,parameters),...
    [parameters.tStart parameters.tFinal], parameters.ros0, options);
[~,time_Dead] = fcn_2; % retrive cell death events

%% The simulated wave plot
% Set polar coordinate (r,phi) for the plotting
r = (0:parameters.nCell-1).*(parameters.h+16); % assume the diameter of a cell is 16 µm
nphi = parameters.nCell-2;
delta_phi = 2*pi/nphi;
phi = (0:parameters.nCell-1).*delta_phi;
[Phi,R] = meshgrid(phi,r);
[x,y] = pol2cart(Phi,R);

% Set colormap for ROS level
load("custom_parula.mat")
CustomColormap(1:2,:) = repmat([0 0 0],2,1); % black color for dead cells

% Visualize the snapshot of the wave at different time points (tp)
figure('Position',[550         450        1100         225])
for tp = 0:30:120

    % Get the ROS level of all cells at a specific time point
    ROS = deval(sol,tp); 
    ROS = reshape(ROS,parameters.nCell,[]);
    
    %  If the cells are dead at the time point, set their ROS to be 0
    if ~isempty(time_Dead)
        ROS(time_Dead(time_Dead(:,2)<tp-parameters.delay,1)) = 0;
    end
    
    % plot ROS wave snapshot at a specific time point
    nexttile
    surf(x, y, ROS,'EdgeColor','none');
    
    % Set format and annotation for the plot    
    axis equal
    set(gca,'BoxStyle','full','Layer','top','LineWidth',1,'view',[0 90],...
        'XLim',1000.*[-1 1],'XTick',[],...
        'YLim',1000.*[-1 1],'YTick',[],'Box','on')
    colormap(CustomColormap)
    caxis([0 4.5])
    title("T = "+tp+" min")
        
end
sgtitle("Erastin = "+parameters.E+" µM");

%% Functions
%% Calculate steady states (SS) of ROS
function SS = steadystates(parameters)
% Solve a symbolic equation with the pre-specified parameter values
% (parameters) to get the steady states (SS)

% Declare the symbolic variables
syms ros k_pfb EC50_pfb n_pfb k_dnfb EC50_dnfb n_dnfb E c_GSH EC50_erastin k_deg k_synth 

% The reaction term of the model
f(ros) = k_pfb*(ros^n_pfb/(EC50_pfb^n_pfb+ros^n_pfb))...
    - k_dnfb*(c_GSH+EC50_erastin/(EC50_erastin+E))*(EC50_dnfb^n_dnfb/(ros^n_dnfb+EC50_dnfb^n_dnfb))*ros...
    - k_deg*ros + k_synth;

% Substitute the values of parameters into the reaction term, f(ros)
f = subs(f,{k_pfb EC50_pfb n_pfb k_dnfb EC50_dnfb n_dnfb E c_GSH EC50_erastin k_deg k_synth},...
    {parameters.k_pfb,parameters.EC50_pfb,parameters.n_pfb,...
    parameters.k_dnfb,parameters.EC50_dnfb,parameters.n_dnfb,...
    parameters.E,parameters.c_GSH,parameters.EC50_erastin,...
    parameters.k_deg,parameters.k_synth});

% Solve f(ros)==0 to get steady states
SS = double(vpasolve(f,ros,[0 Inf]));
end

%% ROS trigger wave model in 2D (polar coordinate)
function dros = ftw2D(t,y,parameters)
% t: current time
% y: ROS level at the current time
% parameters: pre-specified values for parameters

N = parameters.nCell^2; % total cell number in simulation

% Set parameters for individual cells
k_pfb = parameters.k_pfb*ones(N,1);
EC50_pfb = parameters.EC50_pfb*ones(N,1);
n_pfb = parameters.n_pfb*ones(N,1);

k_dnfb = parameters.k_dnfb*ones(N,1);
EC50_dnfb = parameters.EC50_dnfb*ones(N,1);
n_dnfb = parameters.n_dnfb*ones(N,1);

c_GSH = parameters.c_GSH*ones(N,1);
EC50_erastin = parameters.EC50_erastin*ones(N,1);
E = parameters.E*ones(N,1);

k_deg = parameters.k_deg*ones(N,1);
k_synth = parameters.k_synth*ones(N,1);

%% Set parameters for the dead cells
global idx_Dead time_Dead

% Remove undead cells when time jumps backward
if ~isempty(time_Dead)
    idx_Dead(time_Dead(:,2)>t) = [];
    time_Dead(time_Dead(:,2)>t,:) = [];
end

% If there are cells whose ROS level is above the death threshold, record
% it in the idx_Dead
idx_passingTh = find(y-parameters.death_threshold>0);
idx_Dead_new = setdiff(idx_passingTh,idx_Dead);
idx_Dead = [idx_Dead,idx_Dead_new'];

% If there are cells whose ROS level is above the deadth threshold for 30
% min, set k_pfb, k_dnfb, and k_synth to be zero, representing cell death.
if numel(idx_Dead)>=1
    
    % Record the death time
    time_Dead = [time_Dead;idx_Dead_new,t.*ones(length(idx_Dead_new),1)];
    
    % Set the parameters to be zero if cells are dead
    idx = find(time_Dead(:,2)+parameters.delay<t);
    if ~isempty(idx)
        k_pfb(time_Dead(idx,1)) = 0;
        k_dnfb(time_Dead(idx,1)) = 0;
        k_synth(time_Dead(idx,1)) = 0;
    end
end

%% Reshape the ROS level into a matrix form and get the coordinate of the cells (r,phi)
y = reshape(y,parameters.nCell,[]);
nr = size(y,1)-1;
delta_r = parameters.h;
r = ((1:nr)-1).*delta_r;
nphi = size(y,2)-2;
delta_phi = 2*pi/nphi;

% Calculate ROS dynamics for each cell
dros = nan(size(y));
for j = 1:nphi+2
    for i = 1:nr-1
        ind = sub2ind(size(y),i,j);
        
        % ROS dynamics of the cell at r=0
        if i==1
            dros(i,j) = (parameters.D/delta_r^2)*(y(i+1,j)-y(i,j))...
                + (k_pfb(ind)*y(i,j)^n_pfb(ind))/(y(i,j)^n_pfb(ind)+EC50_pfb(ind)^n_pfb(ind))...
                - k_dnfb(ind)*(c_GSH(ind)+EC50_erastin(ind)/(EC50_erastin(ind)+E(ind)))*(EC50_dnfb(ind)^n_dnfb(ind)/(y(i,j)^n_dnfb(ind)+EC50_dnfb(ind)^n_dnfb(ind)))*y(i,j)...
                - k_deg(ind)*y(i,j) + k_synth(ind);

        % ROS dynamics of the cells at phi==0
        elseif j==1
            dros(i,j) = parameters.D*((1/delta_r^2)*(y(i-1,j)-2*y(i,j)+y(i+1,j))...
                + (1/r(i))*(1/(2*delta_r))*(y(i+1,j)-y(i-1,j))...
                + (1/r(i)^2)*(1/delta_phi^2)*(y(i,nphi+1)-2*y(i,j)+y(i,j+1)))...
                + (k_pfb(ind)*y(i,j)^n_pfb(ind))/(y(i,j)^n_pfb(ind)+EC50_pfb(ind)^n_pfb(ind))...
                - k_dnfb(ind)*(c_GSH(ind)+EC50_erastin(ind)/(EC50_erastin(ind)+E(ind)))*(EC50_dnfb(ind)^n_dnfb(ind)/(y(i,j)^n_dnfb(ind)+EC50_dnfb(ind)^n_dnfb(ind)))*y(i,j)...
                - k_deg(ind)*y(i,j) + k_synth(ind);

        % ROS dynamics of the cells at phi==2*pi
        elseif j==nphi+2
            dros(i,j) = parameters.D*((1/delta_r^2)*(y(i-1,j)-2*y(i,j)+y(i+1,j))...
                + (1/r(i))*(1/(2*delta_r))*(y(i+1,j)-y(i-1,j))...
                + (1/r(i)^2)*(1/delta_phi^2)*(y(i,j-1)-2*y(i,j)+y(i,2)))...
                + (k_pfb(ind)*y(i,j)^n_pfb(ind))/(y(i,j)^n_pfb(ind)+EC50_pfb(ind)^n_pfb(ind))...
                - k_dnfb(ind)*(c_GSH(ind)+EC50_erastin(ind)/(EC50_erastin(ind)+E(ind)))*(EC50_dnfb(ind)^n_dnfb(ind)/(y(i,j)^n_dnfb(ind)+EC50_dnfb(ind)^n_dnfb(ind)))*y(i,j)...
                - k_deg(ind)*y(i,j) + k_synth(ind);

        % ROS dynamics of the cells at other positions
        else
            dros(i,j) = parameters.D*((1/delta_r^2)*(y(i-1,j)-2*y(i,j)+y(i+1,j))...
                + (1/r(i))*(1/(2*delta_r))*(y(i+1,j)-y(i-1,j))...
                + (1/r(i)^2)*(1/delta_phi^2)*(y(i,j-1)-2*y(i,j)+y(i,j+1)))...
                + (k_pfb(ind)*y(i,j)^n_pfb(ind))/(y(i,j)^n_pfb(ind)+EC50_pfb(ind)^n_pfb(ind))...
                - k_dnfb(ind)*(c_GSH(ind)+EC50_erastin(ind)/(EC50_erastin(ind)+E(ind)))*(EC50_dnfb(ind)^n_dnfb(ind)/(y(i,j)^n_dnfb(ind)+EC50_dnfb(ind)^n_dnfb(ind)))*y(i,j)...
                - k_deg(ind)*y(i,j) + k_synth(ind);
        end
    end
end
dros(isnan(dros)) = 0;
dros = dros(:); % Convert matrix into a vector form for the usage of ode45
end

%% Record cell death events
function fcn_1(x,y)
global idx_Dead time_Dead
idx_Dead = x;
time_Dead = y;
end

function [x,y] = fcn_2()
global idx_Dead time_Dead
x = idx_Dead;
y = time_Dead;
end

