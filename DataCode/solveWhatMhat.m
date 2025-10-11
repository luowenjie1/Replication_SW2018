function diff = solveWhatMhat(guess, baseline, shocks, loop_shock, n, N, J, parameter)    

    %=============================================================
    % Purpose: Solve the nonlinear system of Eq.(12) + Eq.(13)
    % Inputs:
    %   guess      = Initial guesses (hats of wages and firm numbers)
    %   baseline   = Baseline data (year 1990)
    %   shocks     = Various shocks (Gamma_hat_foreign, Gamma_hat_US, etc.)
    %   loop_shock = Type of shock (1=Foreign, 2=US, 3=Regulation, 4=Expenditure)
    %   n          = Year index
    %   N          = Number of countries
    %   J          = Number of industries
    %   parameter  = Model parameters (alpha, sigma, theta)
    % Output:
    %   diff = Residuals of Eq.(12) + Eq.(13) (passed to fsolve to approach zero)
    %=============================================================


    %% -------- 1. Initialize shock variables (Appendix Step 1) --------
    Gamma_hat = ones(N,N,J); % bilateral productivity shifter
    beta_hat  = ones(N,J);   % expenditure share shock
    NX_hat    = ones(N,1);   % net exports shock (agg.)
    NXs_hat   = ones(N,J);   % net exports shock (sectoral)

    if loop_shock == 1
        Gamma_hat = shocks.Gamma_hat_foreign(:,:,:,n);   % Foreign competitiveness shock
    elseif loop_shock == 2
        Gamma_hat = shocks.Gamma_hat_US(:,:,:,n);        % U.S. competitiveness shock
    elseif loop_shock == 3
        Gamma_hat = shocks.t_hat(:,:,:,n);               % Environmental regulation shock
    elseif loop_shock == 4
        beta_hat  = squeeze(shocks.beta_hat(:,:,:,n));   % Expenditure share shock
    end

    %% -------- 2. Recover wages & firm numbers from guess --------
    % guess = [ w_US ; M_hat(N*J) ]
    w_US = guess(1);

    % Foreign wage determined by numeraire constraint:
    % w_US * (R_US / total R) + w_F * (R_F / total R) = 1
    w_foreign = (1 - w_US .* baseline.R(2, :)) ./ baseline.R(1, :);

    w_hat = [w_foreign; w_US];   % N×1

    % Firm number change (N×J)
    M_hat = reshape(guess(2:end), [N J]);

    %% -------- 3. Eq.(13): Wage condition --------

    w_hat_nnj   = reshape(w_hat, [N 1 1]);       
    w_d_hat_nnj = reshape(w_hat, [1 N 1]);       
    M_hat_nnj   = reshape(M_hat, [N 1 J]);      
    beta_hat_nnj= reshape(beta_hat,[1 N J]);
    NX_d_hat_nnj = reshape(NX_hat, [N 1 1]);

    baseline.R_d_nnj = reshape(baseline.R, [1 N 1]);
    baseline.NX_nj = reshape(baseline.NX, [N 1]);
    baseline.NX_d_nnj = reshape(baseline.NX, [1 N 1]);

    parameter.theta_nnj = reshape(parameter.theta, [1 1 J]);

    diff1 = w_hat_nnj - sum(     (baseline.zeta .* w_hat_nnj.^(-parameter.theta_nnj) .* Gamma_hat) ...
                          ./ (sum(baseline.lambda .* M_hat_nnj .* w_hat_nnj.^(-parameter.theta_nnj) .* Gamma_hat, 1)) ...
                          .* beta_hat_nnj ...
                          .* (    (w_d_hat_nnj .* baseline.R_d_nnj - NX_d_hat_nnj .* baseline.NX_d_nnj) ...
                               ./ (baseline.R_d_nnj - baseline.NX_d_nnj) ...
                              ) ...
                        ,2);

    diff1 = squeeze(diff1);
    diff1 = reshape(diff1, [N*J 1]);

    %% -------- 4. Eq.(12): Market equilibrium condition (resource condition) --------
    % Parameter constants
    cons1 = (parameter.theta - (parameter.sigma-1).*(1-parameter.alpha)) ...
           ./ (parameter.sigma .* parameter.theta);
    cons2 = ((parameter.sigma-1).*(parameter.theta-parameter.alpha+1)) ...
           ./ (parameter.sigma .* parameter.theta);

    cons1 = reshape(cons1, [1 J]);
    cons2 = reshape(cons2, [1 J]);

    % ψ:  (Footnote 16)
    psi =  (1 - sum(cons1 .* baseline.beta, 2)) ...
        ./ (1 - sum(cons1 .* beta_hat .* baseline.beta, 2));
    psi = reshape(psi, [N 1]);

    % η:  (Footnote 16)
    eta_0 = sum( (-cons1+1) .* baseline.beta .* baseline.NX_nj ...
                 - baseline.NXs .* cons2 , 2);

    % η_post: 
    eta_0_post = sum( (-cons1+1).*beta_hat.*baseline.beta.*NX_hat.*baseline.NX_nj ...
                      - NXs_hat.*baseline.NXs.*cons2, 2);
    eta_0_post = eta_0_post ./ w_hat;


    eta_0      = reshape(eta_0, [N 1]);
    eta_0_post = reshape(eta_0_post, [N 1]);

    % Eq.(12): resource condition
    diff2 = 1 - psi .* ( sum(M_hat .* baseline.Rs .* cons2, 2) + eta_0_post ) ...
                    ./ ( sum(baseline.Rs .* cons2, 2) + eta_0 );
    diff2 = diff2(2:end);   % drop one equation (numeraire country)

    %% -------- 5. Combine residuals of Eq.(12)+(13) --------
    diff = [diff2; diff1];

end
