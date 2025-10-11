%======================================================================
% Date:        2025-10-10
% Author:      Wenjie Luo
%
% Purpose:
%   This script solves the nonlinear system consisting of Eq. (12) and 
%   Eq. (13) in the model. The system determines the equilibrium wage 
%   and the number of firms (in hat values) under various counterfactual
%   scenarios. Numerical methods (fsolve) are used to find the solution 
%   based on the baseline calibration and specified shocks.
%
% Reference:
%   - Shapiro, J.S. and Walker, R. (2018). “Why Is Pollution from U.S. 
%     Manufacturing Declining? The Roles of Environmental Regulation, 
%     Productivity, and Trade.” *American Economic Review*, 108(12), 
%     3814–3854.
%
% Notes:
%   - Input data "rawFile.mat" comes from replication file in SW2018
%   - CHANGE "pollution_index" \in [1, 7] in Line 39 for different pollutants
%======================================================================

clear; 
clc;

% Try to set the working directory automatically
if isempty(mfilename('fullpath'))
    % Fall back to current folder if script is run interactively
    currentPath = pwd;
else
    currentPath = fileparts(mfilename('fullpath'));
    cd(currentPath);
end

%% ===================== Basic Settings =====================

% Select the index of the pollutant (1=CO, 2=NOx, 3=PM10, …)
pollutants = {'co','nox','pm10','pm25','so2','voc','co2'};
pollution_index = 1;

% for pollutant = 1:length(pollutants)
% pollution_index = pollutant;
% disp(pollutants(pollution_index));

% load data
load rawFile.mat

% Data structure: xbilat(o,d,j,t)
% o: origin countries (2: foreign, US)
% d: destination countries (2: foreign, US)
% j: industries (17 industries)
% t: years (19 years, 1990–2008)
%
% Original data:
% co2_poll(j,t) - CO2 emissions (6 observed years)
% poll(j,t,p) - emissions of 7 pollutants (6 observed years)
% vship(o,j,t) - production by country and industry (equivalent to Rs)
% xbilat(o,d,j,t) - bilateral trade flows

%% ===================== Parameter Settings =====================
N  = 2;       % Number of countries (1=foreign, 2=US)
us = 2;       % US index
J  = 17;      % Number of industries
Y  = 19;      % Number of years
yvec = 1990:2008;

% run plot_pollution_trend.m

% parameter: model parameters
% alpha: pollution elasticity: Table 2, Column (2)
parameter.alpha = [.0040; .0022; .0103; .0223; .0212; .0205; .0048; .0303; .0557; .0019; .0015; .0023; .0005; .0014; .0016; .0019; .0047];

% inputshare: wL/R: Table 2, Column (3)
inputshare = [.74;.79;.83;.79;.88;.70;.78;.73;.85;.79;.76;.81;.79;.65; .82;.74;.73];

% sigma: elasticity of substitution <- Footnote 8, App p19
parameter.sigma = (1-parameter.alpha)./((1-parameter.alpha)-inputshare);  

% theta: Pareto shape parameter: Table 2, Column (5)
parameter.theta = [4.81; 5.38; 8.30; 4.29; 17.52; 4.13; 5.02; 3.39; 9.72; 5.60; 4.30; 5.07; 4.13; 2.09; 5.29; 3.27; 4.77];      

% pm.dirty: dirty industries; pm.dirty: clean industries
% parameter.dirty = [4 5 6 8 9];
% parameter.clean = [1 2 3 7 10 11 12 13 14 15 16 17];

%% ===================== Data Processing =====================
% Total import expenditure and total revenue
Es = squeeze(sum(xbilat,1));     % destination-sector-year (∑origin)
Rs = squeeze(sum(xbilat,2));     % origin-sector-year (∑destination)
NXs = Rs - Es;                   % Net exports by industry

E = squeeze(sum(Es,2));          % destination-year
R = squeeze(sum(Rs,2));          % origin-year
NX = R - E;                      % Net exports by country

% β: CES expenditure share
beta   = Es ./ reshape(E, [N 1 Y]);             % d×j×t
lambda = xbilat ./ reshape(Es, [1 N J Y]);      % o×d×j×t, expenditure share
zeta   = xbilat ./ reshape(Rs, [N 1 J Y]);      % o×d×j×t, revenue share

% extrapolate CO2 emissions
%   [1991,1994,1998,2002,2006,2010] -> (1990:2008)
for j = 1:J
    Z_co2(j,:) = interp1([1991,1994,1998,2002,2006,2010]',co2_poll(j,:)',(1990:2008)','linear','extrap')';
    for n = 1:length(pollutants)
        Z(j,:,n) = interp1q([1990,1996,1999,2002,2005,2008]',poll(j,:,n)',(1990:2008)')';
    end
end
Z(:,:,7) = Z_co2;

%% ===================== Baseline Data (t=1, 1990) =====================
baseline.xbilat = xbilat(:,:,:,1);
baseline.lambda = lambda(:,:,:,1);
baseline.zeta   = zeta(:,:,:,1);
baseline.Rs     = Rs(:,:,1);        % origin×sector
baseline.R      = R(:,1);           % origin
baseline.NXs    = NXs(:,:,1);       % origin×sector
baseline.NX     = NX(:,1);          % origin
baseline.Z      = Z(:,1,:);         % sector×pollutant
baseline.wage   = baseline.R;       % Baseline wage ∝ income
baseline.beta   = beta(:,:,1);      % destination×sector

%% ===================== Hat Variables =====================
lambda_hat = lambda ./ baseline.lambda;
Z_hat      = Z ./ baseline.Z;
Rs_hat     = Rs ./ baseline.Rs;
R_hat      = R ./ baseline.R;
NXs_hat    = NXs ./ baseline.NXs;
NX_hat     = NX ./ baseline.NX;
w_hat      = R ./ baseline.wage;                   % Equation (13), wage hats
M_hat      = Rs_hat ./ reshape(w_hat, [N 1 Y]);    % Equation (11), mass of firms hats

%% ===================== Shocks =====================
beta_hat = beta ./ baseline.beta;

% Pollution tax (Eq. 15)
t_hat = reshape(M_hat, [N 1 J Y]) ...
      .* reshape(w_hat, [N 1 1 Y]) ...
      ./ reshape(Z_hat(:,:,pollution_index), [1 1 J Y]);
t_hat = squeeze(t_hat);

implicit_tax(pollution_index,:) = squeeze(mean(t_hat(us,:,:),2));

% Construct 4D arrays (o×d×j×t) to be substituted into Eq.19 / Eq.21
M_hat_4D     = reshape(M_hat, [N 1 J Y]);
w_hat_4D     = reshape(w_hat, [N 1 1 Y]);
w_d_hat_4D   = reshape(w_hat, [1 N 1 Y]);
beta_hat_4D  = reshape(beta_hat, [1 N J Y]);
R_hat_4D     = reshape(R_hat, [N 1 1 Y]);
R_d_hat_4D = reshape(R_hat, [1 N 1 Y]);
NX_hat_4D    = reshape(NX_hat,[N 1 1 Y]);
NX_d_hat_4D = reshape(NX_hat, [1 N 1 Y]);
t_hat_4D = reshape(t_hat, [N 1 J Y]);

baseline.R_d_4D = reshape(baseline.R, [1 N 1 1]);
baseline.NX_d_4D = reshape(baseline.NX,[1 N 1 1]);

parameter.theta_4D = reshape(parameter.theta,[1 1 J 1]);
parameter.sigma_4D = reshape(parameter.sigma,[1 1 J 1]);
parameter.alpha_4D = reshape(parameter.alpha,[1 1 J 1]);


%% ===================== Gamma (Eq.19, Eq.21) =====================
Gamma_hat_star = lambda_hat ...
               ./ (M_hat_4D .* w_hat_4D.^(-parameter.theta_4D)) ...
               .* ( (beta_hat_4D ./ w_d_hat_4D) ...
                  .* ((R_d_hat_4D.*baseline.R_d_4D - NX_d_hat_4D.*baseline.NX_d_4D) ...
                    ./ (baseline.R_d_4D - baseline.NX_d_4D)) ...
                  ).^( 1 - parameter.theta_4D ./ ((parameter.sigma_4D - 1).*(1 - parameter.alpha_4D)) );
% Foreign only
Gamma_hat_foreign = [Gamma_hat_star(1,:,:,:); ones(1,2,J,Y)];
% US only, with environmental tax shock (Eq.21)
Gamma_hat_US = Gamma_hat_star ./ t_hat_4D.^( -(parameter.alpha_4D.*parameter.theta_4D)./(1-parameter.alpha_4D) );
Gamma_hat_US = [ones(1,2,J,Y); Gamma_hat_US(2,:,:,:)];

% Save
shocks.Gamma_hat_foreign = Gamma_hat_foreign;
shocks.Gamma_hat_US = Gamma_hat_US;
shocks.beta_hat = reshape(beta_hat,[1 N J Y]);

t_hat_shcok = t_hat_4D.^( - (parameter.alpha_4D .* parameter.theta_4D) ...
                   ./ (1 - parameter.alpha_4D) ...
                   );
t_hat_shcok = squeeze(t_hat_shcok);
shocks.t_hat = reshape([ones(1, 17, 19); t_hat_shcok(2, :, :)], [N 1 17 19]);


%% ===================== Step 1: Construct Initial Guess =====================
% Corresponds to Appendix p.22 Step 1
% Initial guess: current wage change (w_hat) and firm number change (M_hat)
% Dimension: [ (N−1 wages) + (N*J firms) ] × 1 vector
for n = 1:Y
    guess(:, n) = [ w_hat(2, n); ...                         % US wage change (foreign wage as numeraire)
                    reshape(M_hat(:,:,n), [N*J, 1]) ];       % Firm number change for each country×industry
end

%% ===================== Step 2: Solve the Nonlinear System Using fsolve =====================
% Corresponds to Appendix p.22 Step 2
% Solve the system of equations (Eq.12, Eq.13), unknowns = w (N−1) + M (N*J)
% Number of equations = number of unknowns → system solvable
options = optimset('Display','off', ...        
                   'MaxFunEvals',60000, ...
                   'MaxIter',4500, ...
                   'TolFun',1e-14, ...
                   'TolX',1e-14, ...
                   'Algorithm','trust-region-dogleg');

for n = 1:Y
    for loop_shock = 1:4   % Four types of shocks (Appendix p.22 Step 3)
        initial_guess = squeeze(guess(:, n));

        % Call fsolve to solve Eq.(12)+(13)
        [solver, fval, flag] = fsolve(@(g) solveWhatMhat(g, baseline, shocks, loop_shock, n, N, J, parameter), ...
                                      initial_guess, options);

        % ================== Step 3: Update Variables (Eq.12, Eq.13) ==================
        % Results from fsolve include:
        % solver(1) = w_US (wage change)
        % solver(2:end) = M_hat (firm number change, reshaped to N×J)

        % Recover wage changes (foreign wage determined by numeraire)
        w_US = solver(1);       
        w_foreign = (1 - w_US .* (baseline.R(2) ./ sum(baseline.R))) ...
                  ./ (baseline.R(1) ./ sum(baseline.R));
        % w_us * US's GDP(or wage) share + w_foreign * Foreign's GDP(or wage) share = 1
        w_hat = [w_foreign; w_US];  % N×1 vector

        % Firm number change (N×J matrix)
        M_hat = reshape(solver(N:end), [N J]);

        % Expand wages to N×J for emission calculation
        w_hat_nj = reshape(w_hat,[N 1]);

        % ================== Counterfactual Emissions (Eq.15, Eq.19, Eq.21) ==================
        if loop_shock == 3
            % Regulation shock
            Z_hat_counterfactual(:,:,n,loop_shock) = (w_hat_nj .* M_hat) ./ squeeze(t_hat(:,:,n));
        else
            % Other shocks, pollution tax fixed at baseline
            Z_hat_counterfactual(:,:,n,loop_shock) = (w_hat_nj .* M_hat) ./ squeeze(t_hat(:,:,1));
        end

        %================== EXTENSION: Counterfactual Revenue ==================
        R_s_counterfactual(:,:,n,loop_shock) = (w_hat_nj .* M_hat);
    end
end

%% ===================== Calculate Aggregate Emissions (Weighted Average, Eq.22) =====================
% Weight = baseline pollution level (Z0)
Z0_sum = sum(baseline.Z(:,:,pollution_index),1);  % Sectoral baseline pollution ∑origin
Z0_sum = reshape(Z0_sum, [1 1 1]);

% Counterfactual US emission path (weighted sum, converted to %)
Z_hat_counterfactual_us = squeeze(Z_hat_counterfactual(us,:,:,:));
Z_hat_counterfactual_us = sum(Z_hat_counterfactual_us .* baseline.Z(:,:,pollution_index),1) ./ Z0_sum;
Z_hat_counterfactual_us = squeeze(Z_hat_counterfactual_us) .* 100;

% Actual data (observed emission trajectory)
Z_hat_real = squeeze(sum(Z(:,:,pollution_index),1))';
Z_hat_real = Z_hat_real ./ Z_hat_real(1) .* 100;

% Combine actual and counterfactual results
Z_hat_counterfactual_us = [Z_hat_real'; Z_hat_counterfactual_us']';

%% ===================== EXTENSION: Calculate Aggregate Revenue =====================
% Counterfactual US emission path (weighted sum, converted to %)
R_s_hat_counterfactual_us = squeeze(R_s_counterfactual(us,:,:,:));
R_hat_counterfactual_us = squeeze(mean(R_s_hat_counterfactual_us, 1));
R_hat_counterfactual_us = R_hat_counterfactual_us .* baseline.R(us);
R_hat_counterfactual_us = R_hat_counterfactual_us ./ R_hat_counterfactual_us(1) .* 100;

% Actual data (observed emission trajectory)
R_hat_real = R(us,:) ./ R(us,1) .* 100;

% Combine actual and counterfactual results
R_hat_real = [R_hat_real; R_hat_counterfactual_us']';

%% ===================== Plotting =====================
% Create a folder called "figure" (if it doesn't exist)
figDir = fullfile(currentPath, 'figure');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

nei_yrs = [1990 1996 1999 2002 2005 2008];
nei_yrsIndex = [1 7 10 13 16 19];

% Pollution trend
f = figure('Visible','on');
figHandles(pollution_index) = f; % store handle for later combination
clf;

plot(1990:2008,Z_hat_counterfactual_us(:,1),'-','LineWidth',2.5); hold on;    

plot(nei_yrs,Z_hat_counterfactual_us(nei_yrsIndex,2),'rp','LineWidth',1.5); 
plot(nei_yrs,Z_hat_counterfactual_us(nei_yrsIndex,3),'ro','LineWidth',1.5); 
plot(nei_yrs,Z_hat_counterfactual_us(nei_yrsIndex,4),'rs','LineWidth',1.5); 
plot(nei_yrs,Z_hat_counterfactual_us(nei_yrsIndex,5),'rv','LineWidth',1.5); 

plot(1990:2008,Z_hat_counterfactual_us(:,2),'--r','LineWidth',1.5); 
plot(1990:2008,Z_hat_counterfactual_us(:,3),'--r','LineWidth',1.5); 
plot(1990:2008,Z_hat_counterfactual_us(:,4),'--r','LineWidth',1.5); 
plot(1990:2008,Z_hat_counterfactual_us(:,5),'--r','LineWidth',1.5); 

set(gca,'FontSize',13,'YTick',[0 30 60 90 120 150],'XTick',[1990 1995 2000 2005 2010]); 
axis([1990 2008 0 150]); 
xlabel('Year'); ylabel('1990=100'); box off;
legend('Actual Data (All Shocks)', ...
           'Foreign Competitiveness Shocks Only', ...
           'U.S. Competitiveness Shocks Only', ...
           'U.S. Regulation Shocks Only', ...
           'U.S. Expenditure Share Shocks Only', ...
           'Location','Southwest');
title(upper(pollutants{pollution_index}));
saveas(f, fullfile(figDir, [pollutants{pollution_index}, '.png']));
close all;

% EXTENSION: revenue trend
if pollution_index == 1
    f = figure('Visible','on');
    clf;
    
    plot(1990:2008,R_hat_real(:,1),'-','LineWidth',2.5); hold on;    
    
    plot(nei_yrs,R_hat_real(nei_yrsIndex,2),'rp','LineWidth',1.5); 
    plot(nei_yrs,R_hat_real(nei_yrsIndex,3),'ro','LineWidth',1.5); 
    plot(nei_yrs,R_hat_real(nei_yrsIndex,4),'rs','LineWidth',1.5); 
    plot(nei_yrs,R_hat_real(nei_yrsIndex,5),'rv','LineWidth',1.5); 
    
    plot(1990:2008,R_hat_real(:,2),'--r','LineWidth',1.5); 
    plot(1990:2008,R_hat_real(:,3),'--r','LineWidth',1.5); 
    plot(1990:2008,R_hat_real(:,4),'--r','LineWidth',1.5); 
    plot(1990:2008,R_hat_real(:,5),'--r','LineWidth',1.5); 
    
    set(gca,'FontSize',13,'YTick',[0 30 60 90 120 150],'XTick',[1990 1995 2000 2005 2010]); 
    axis([1990 2008 0 150]); 
    xlabel('Year'); ylabel('1990=100'); box off;
    legend('Actual Data (All Shocks)', ...
               'Foreign Competitiveness Shocks Only', ...
               'U.S. Competitiveness Shocks Only', ...
               'U.S. Regulation Shocks Only', ...
               'U.S. Expenditure Share Shocks Only', ...
               'Location','Southwest');
    title('US revenue share of the world');
    saveas(f, fullfile(figDir, 'revenue.png'));
    close all;
end

end

%% ===================== Combine All Pollutant Plots =====================
% change `combine_plots` to one once finishing the seven pollutants
combine_plots = 0;
if combine_plots == 1
    run combine_plots.m
end