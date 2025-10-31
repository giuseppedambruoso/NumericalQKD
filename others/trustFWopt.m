% qstar_vs_alpha_opt.m
% Max continuo del rate SOLO Step 2 su q ∈ [0, 0.1] con fminbnd.
% r(q) = f_lb_step2(beta,p,q) - h2(p+q-2pq)  (clip a 0 solo alla fine)
%
% Output:
%   - Curve q*(alpha) per vari p
%   - Heatmap di q*(alpha,p)
%   - Heatmap del rate massimo r*(alpha,p)
%   - ARRAY 3D DATA(numAlpha, numP, 5) con i valori a q*:
%       1: rate_clipped, 2: rate_unc, 3: f_lb, 4: leak, 5: q_star
% Richiede: bb84_fbeta_min.m (Step 1+2)

clear all; clc;

% ---- set esperimento ----
alpha_list = [1.0001, 1.025, 1.05, 1.075, 1.10, 1.125, 1.15, 1.175, 1.20, ...
    1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3];
beta_list  = 1 ./ alpha_list; %#ok<NASGU>
p_list     = [0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10, 0.105, 0.11, 0.115, 0.12];

q_lo = 0.0; 
q_hi = 0.5;                                  % dominio di ottimizzazione su q
opts_fmin = optimset('TolX', 1e-3, 'Display', 'off');

% Opzioni per Frank–Wolfe interno (Step 1)
optsFW = struct('maxIter',40,'tol_gap',1e-5,'eps_reg',1e-10,'verbose',false);

% ---- contenitori ----
numA = numel(alpha_list);
numP = numel(p_list);

Qstar = nan(numA, numP);    % q* per (alpha,p)
Rstar = nan(numA, numP);    % r*(alpha,p) = max_q r(q)  (clippato a 0)
Runc  = nan(numA, numP);    % r_unc*(alpha,p) = f_lb - h2 senza clip
R0    = nan(numA, numP);    % r(q=0) clippato
R0unc = nan(numA, numP);    % r(q=0) non clippato
F_LB  = nan(numA, numP);    % f_lb_step2(beta,p,q*) a q*
LEAK  = nan(numA, numP);    % leak h2(p+q*-2pq*) a q*

% -------- Array 3D con tutti i valori a q* --------
% Slices:
% 1: rate_clipped, 2: rate_unc, 3: f_lb, 4: leak, 5: q_star
DATA  = nan(numA, numP, 5);

for ip = 1:numP
    p = p_list(ip);
    fprintf('\n=== p = %.3f ===\n', p);
    for ia = 1:numA
        tic;
        alpha = alpha_list(ia);
        beta  = 1 / alpha;

        % Definisci il rate senza clip come funzione di q per l'ottimizzatore
        rate_unc = @(q) local_rate_unc(beta, p, q, optsFW);  % f_lb_step2 - h2(...)

        % --- valore a q=0 (baseline) ---
        rq_0_unc = rate_unc(q_lo);
        r0       = max(rq_0_unc, 0);
        R0unc(ia, ip) = rq_0_unc;
        R0(ia, ip)    = r0;

        % Massimizza con fminbnd minimizzando il negativo
        obj = @(q) -rate_unc(q);
        [q_cont, neg_val, ~] = fminbnd(obj, q_lo, q_hi, opts_fmin);
        rq_cont = -neg_val;

        % Confronta anche gli estremi del dominio (robustezza)
        rq_0 = rate_unc(q_lo);
        rq_1 = rate_unc(q_hi);

        % Scegli il migliore tra {q_cont, q_lo, q_hi}; tie-break: q più piccolo
        r_candidates = [rq_cont, rq_0, rq_1];
        q_candidates = [q_cont, q_lo, q_hi];
        [r_best_unc, idx_best] = max(r_candidates);
        q_best = q_candidates(idx_best);

        % Valori esatti a q_best
        leak_best = h2(p + q_best - 2*p*q_best);
        f_lb_best = r_best_unc + leak_best; 
        r_best    = max(r_best_unc, 0);

        % Salva in matrici 2D
        Qstar(ia, ip) = q_best;
        Rstar(ia, ip) = r_best;
        Runc(ia, ip)  = r_best_unc;
        %F_LB(ia, ip)  = f_lb_best;
        %LEAK(ia, ip)  = leak_best;

        % % Salva nel CUBO 3D
        % DATA(ia, ip, 1) = r_best;      % rate clippato
        % DATA(ia, ip, 2) = r_best_unc;  % rate non clippato
        % DATA(ia, ip, 3) = f_lb_best;   % f_lb
        % DATA(ia, ip, 4) = leak_best;   % leak
        % DATA(ia, ip, 5) = q_best;      % q*

        fprintf(' alpha=%.4f (beta=%.6f) -> q* = %.4f,  rate* = %.6f (unc=%.6f)\n', ...
                alpha, 1/alpha, q_best, r_best, r_best_unc);
        toc
    end
end


%%

%% ======= Δ r_alpha = max_p | r(p,q=0) - r(p,q*) |  e plot vs alpha =======
Delta = zeros(numA,1);
p_at_max = nan(numA,1);

for ia = 1:numA
    diffs = abs(R0(ia,:) - Rstar(ia,:));   % r(q=0) vs r(q*)
    [Delta(ia), idx] = max(diffs);
    p_at_max(ia) = p_list(idx);
end

% Plot Δ r_alpha vs alpha
figure; hold on; grid on; box on;
plot(alpha_list, Delta, '-o', 'LineWidth', 2, 'MarkerSize', 7);
xlabel('$\alpha$','FontSize',26,'Interpreter','latex');
ylabel('$\Delta r_\alpha$ ','FontSize',26,'FontName','Computer Modern','Interpreter','latex');
title('Massima variazione del rate con trusted noise vs \alpha');
% (opzionale) annota il p che massimizza
for ia = 1:numA
    text(alpha_list(ia), Delta(ia), sprintf(' p=%.3f', p_at_max(ia)), ...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
end

% Log a console
fprintf('\nΔr_alpha e p che la massimizza:\n');
for ia = 1:numA
    fprintf(' alpha=%.2f: Δr = %.6f  (p = %.3f)\n', alpha_list(ia), Delta(ia), p_at_max(ia));
end



%%
% ---- plot q*(alpha) per ciascun p (curve) ----
figure; hold on; grid on; box on;
for ip = 1:numP
    plot(alpha_list, Qstar(:,ip), '-o', 'LineWidth', 2, 'MarkerSize', 6);
end
xlabel('\alpha'); ylabel('q^\star(\alpha)');
title('Trusted-noise ottimale q^\star vs \alpha (rate Step 2, max continuo su q)');
legend(arrayfun(@(pp) sprintf('p=%.2f', pp), p_list, 'UniformOutput', false), ...
       'Location','northwest');

%%

% ---- HEATMAP di q*(alpha,p) ----
figure;
imagesc(alpha_list, p_list, Qstar');   % righe=p, colonne=alpha
set(gca, 'YDir', 'normal');
colorbar; caxis([q_lo q_hi+eps]);
xlabel('\alpha'); ylabel('p');
title('Heatmap di q^\star(\alpha,p)  (rate Step 2, max continuo su q)');
colormap parula;

% Etichette numeriche (opzionali)
hold on;
for ia = 1:numA
    for ip = 1:numP
        text(alpha_list(ia), p_list(ip), sprintf('%.2f', Qstar(ia,ip)), ...
            'HorizontalAlignment','center','Color','w','FontSize',8,'FontWeight','bold');
    end
end
hold off;

% ---- (facoltativo) heatmap del rate massimo r*(alpha,p) ----
figure;
imagesc(alpha_list, p_list, Rstar');   % righe=p, colonne=alpha
set(gca, 'YDir', 'normal');
colorbar;
xlabel('\alpha'); ylabel('p');
title('Heatmap del rate massimo r^*(\alpha,p)  (Step 2, max continuo su q)');
colormap turbo;
%%
% ---- grafici separati per ogni alpha: r(q=0) e r(q*) in funzione di p ----
for ia = 1:numA
    figure; hold on; grid on; box on;
    plot(p_list, R0(ia,:),   '--o', 'LineWidth', 2, 'MarkerSize', 6);
    plot(p_list, Rstar(ia,:), '-s', 'LineWidth', 2, 'MarkerSize', 6);
    xlabel('p'); ylabel('rate');
    ylim([0.0 0.3]);
    title(sprintf('Rate vs p (\\alpha=%.2f, \\beta=%.4f) — Step 2', alpha_list(ia), 1/alpha_list(ia)));
    legend('q = 0', 'q = q^\ast', 'Location','best');
end

%%

surf(alpha_list,p_list,Qstar')

%%

% ======================= funzioni locali ========================

function r_unc = local_rate_unc(beta, p, q, optsFW)
    % Valuta il rate senza clip: r_unc(q) = f_lb_step2(beta,p,q) - h2(p+q-2pq)
    [~, ~, f_lb_step2, ~] = bb84_fbeta_min(beta, p, q, optsFW);
    r_unc = f_lb_step2 - h2(p + q - 2*p*q);
end

function H = h2(x)
    x = min(max(real(x), 1e-15), 1-1e-15);
    H = -x.*log2(x) - (1-x).*log2(1-x);
end



function [f_star, rho_star, f_lb_step2, sigma_step2] = bb84_fbeta_min(beta, p, q, opts)
% Minimizza f_beta (Step 1) e calcola il lower bound lineare (Step 2).
% OUT:
%   f_star       = min f_beta trovato (Step 1)
%   rho_star     = stato ottimo trovato
%   f_lb_step2   = lower bound rigoroso via linearizzazione in rho_star (Step 2)
%   sigma_step2  = argmin della linearizzazione (stato per Step 2)

    if nargin<4, opts=struct; end
    if ~isfield(opts,'maxIter'),  opts.maxIter  = 60;     end
    if ~isfield(opts,'tol_gap'),  opts.tol_gap  = 1e-7;   end
    if ~isfield(opts,'eps_reg'),  opts.eps_reg  = 1e-10;  end
    if ~isfield(opts,'verbose'),  opts.verbose  = false;  end

    % --- basi e proiettori ---
    Id2 = eye(2);
    ket0 = [1;0]; ket1=[0;1];
    plus = (ket0+ket1)/sqrt(2); minus=(ket0-ket1)/sqrt(2);
    Z0 = ket0*ket0'; Z1 = ket1*ket1';
    Pp = plus*plus'; Pm = minus*minus';
    PiZerr = kron(Z0,Z1) + kron(Z1,Z0);
    PiXerr = kron(Pp,Pm) + kron(Pm,Pp);

    % --- mappa G con trusted noise q ---
    sqrtL0 = diag([sqrt(1-q), sqrt(q)]);
    sqrtL1 = diag([sqrt(q),   sqrt(1-q)]);
    K0 = kron(Id2, sqrtL0);
    K1 = kron(Id2, sqrtL1);

    % --- Step 1: punto iniziale fattibile ---
    rho = find_feasible_rho(PiZerr, PiXerr, p);

    % --- Step 1: Frank–Wolfe ---
    f_prev = Inf; f_curr = Inf;
    for it = 1:opts.maxIter
        [fval, grad] = objective_and_gradient(rho, beta, K0, K1, opts.eps_reg);
        sigma = linear_oracle(grad, PiZerr, PiXerr, p);           % direzione FW (SDP lineare)
        gap   = real(trace((rho - sigma)'*grad));                 % dual gap
        lambda = 1.0; f_best = fval; rho_best = rho;              % line search
        while true
            cand = (1-lambda)*rho + lambda*sigma;
            f_new = objective_only(cand, beta, K0, K1, opts.eps_reg);
            if f_new <= f_best - 1e-10 || lambda < 1/1024
                f_best = f_new; rho_best = cand; break;
            end
            lambda = lambda/2;
        end
        rho = rho_best; f_curr = f_best;
        if opts.verbose
            fprintf('[%02d] f=%.12f, gap=%.3e, λ=%.3f\n', it, f_curr, gap, lambda);
        end
        if gap < opts.tol_gap || abs(f_prev - f_curr) < 1e-9, break; end
        f_prev = f_curr;
    end

    f_star   = f_curr;      % valore Step 1
    rho_star = rho;

    % --- Step 2: lower bound tramite linearizzazione in rho_star ---
    [f_at_rho, grad_at_rho] = objective_and_gradient(rho_star, beta, K0, K1, opts.eps_reg);
    sigma_step2 = linear_oracle(grad_at_rho, PiZerr, PiXerr, p); % risolve min <∇f(ρ*), σ>
    % L(σ) = f(ρ*) + <∇f(ρ*), σ-ρ*>,  con ∇f(ρ*) = f*I + G†[...];  Tr(σ)=Tr(ρ*)=1 ⇒ il termine f*I si cancella
    f_lb_step2 = f_at_rho + real(trace(grad_at_rho*(sigma_step2 - rho_star)));

end

% ======================== Helper (identici a prima) ========================
function rho = find_feasible_rho(PiZerr, PiXerr, p)
    cvx_begin quiet sdp
        variable rho(4,4) hermitian semidefinite
        minimize( 0 )
        subject to
            trace(rho) == 1;
            for a=1:2
               for ap=1:2
                   ind = @(b) ( (a-1)*2 + b );
                   indp= @(bp)((ap-1)*2 + bp);
                   sum_b = 0;
                   for b=1:2
                       sum_b = sum_b + rho(ind(b), indp(b));
                   end
                   sum_b == (a==ap)*0.5;
               end
            end
            real(trace(PiZerr*rho)) == p;
            real(trace(PiXerr*rho)) == p;
    cvx_end
    if ~contains(cvx_status,'Solved')
        error('Feasibility SDP failed: %s', cvx_status);
    end
end

function sigma = linear_oracle(grad, PiZerr, PiXerr, p)
    cvx_begin quiet sdp
        variable sigma(4,4) hermitian semidefinite
        minimize( real(trace(grad*sigma)) )
        subject to
            trace(sigma) == 1;
            for a=1:2
               for ap=1:2
                   ind = @(b) ( (a-1)*2 + b );
                   indp= @(bp)((ap-1)*2 + bp);
                   sum_b = 0;
                   for b=1:2
                       sum_b = sum_b + sigma(ind(b), indp(b));
                   end
                   sum_b == (a==ap)*0.5;
               end
            end
            real(trace(PiZerr*sigma)) == p;
            real(trace(PiXerr*sigma)) == p;
    cvx_end
    if ~contains(cvx_status,'Solved')
        warning('Linear oracle SDP status: %s', cvx_status);
    end
end

function [f, grad] = objective_and_gradient(rho, beta, K0, K1, eps_reg)
    [Grho, Sigma] = build_G_and_Zpinch(rho, K0, K1);
    [f, DBgrad]   = renyi_div_and_grad(Grho, Sigma, beta, eps_reg);
    Gdag = @(X) G_dagger_to_AB(X, K0, K1);
    grad = f*eye(4) + Gdag(DBgrad);
end

function f = objective_only(rho, beta, K0, K1, eps_reg)
    [Grho, Sigma] = build_G_and_Zpinch(rho, K0, K1);
    f = renyi_div_only(Grho, Sigma, beta, eps_reg);
end

function [Grho, Sigma] = build_G_and_Zpinch(rho, K0, K1)
    Grho = zeros(8,8);
    K = {K0, K1};
    for j=1:2
        for k=1:2
            bj = (j-1)*4 + (1:4);
            bk = (k-1)*4 + (1:4);
            Grho(bj, bk) = K{j}*rho*K{k};
        end
    end
    Sigma = zeros(8,8);
    Sigma(1:4,1:4) = Grho(1:4,1:4);
    Sigma(5:8,5:8) = Grho(5:8,5:8);
end

function X = G_dagger_to_AB(X8, K0, K1)
    X = zeros(4,4);
    K = {K0, K1};
    for j=1:2
      for k=1:2
          bj = (j-1)*4 + (1:4);
          bk = (k-1)*4 + (1:4);
          X = X + K{j}' * X8(bj,bk) * K{k};
      end
    end
end

function [Dval, GradForward] = renyi_div_and_grad(R, S, beta, eps_reg)
    mu = (1 - beta)/(2*beta);
    R = (R+R')/2; S = (S+S')/2; S = S + eps_reg*eye(size(S));
    [Us, Ds] = eig(S);
    ds = max(real(diag(Ds)), eps_reg);
    S_mu = Us*diag(ds.^mu)*Us';
    Xi = S_mu * R * S_mu; Xi = (Xi+Xi')/2;
    [Ux, Dx] = eig(Xi);
    dx = max(real(diag(Dx)), eps_reg);
    Q = sum(dx.^beta);
    Dval = (1/(beta-1))*log2( Q / max(real(trace(R)), eps_reg) );
    Xi_bm1 = Ux*diag(dx.^(beta-1))*Ux';
    chi2 = beta * (S_mu * Xi_bm1 * S_mu);
    Tmu = @(A) spectral_kernel(A, Us, ds, mu);
    chi1 = beta * Zpinch( Tmu( R * S_mu * Xi_bm1 ) );
    chi3 = beta * Zpinch( Tmu( Xi_bm1 * S_mu * R ) );
    GradForward = ( (chi1+chi2+chi3)/Q - eye(size(R))/max(real(trace(R)),eps_reg) ) / (beta-1);
end

function Dval = renyi_div_only(R, S, beta, eps_reg)
    mu = (1 - beta)/(2*beta);
    R = (R+R')/2; S = (S+S')/2; S = S + eps_reg*eye(size(S));
    [Us, Ds] = eig(S);
    ds = max(real(diag(Ds)), eps_reg);
    S_mu = Us*diag(ds.^mu)*Us';
    Xi = S_mu * R * S_mu; Xi = (Xi+Xi')/2;
    [Ux, Dx] = eig(Xi);
    dx = max(real(diag(Dx)), eps_reg);
    Dval = (1/(beta-1))*log2( sum(dx.^beta) / max(real(trace(R)), eps_reg) );
end

function Y = spectral_kernel(A, U, d, mu)
    tol = 1e-14;
    Atil = U' * A * U;
    n = length(d); M = zeros(n,n);
    for i=1:n
        for j=1:n
            if i==j
                M(i,j) = mu * max(d(i),tol)^(mu-1);
            else
                denom = d(i) - d(j);
                if abs(denom) < tol
                    M(i,j) = mu * max(d(i),tol)^(mu-1);
                else
                    M(i,j) = ( max(d(i),tol)^mu - max(d(j),tol)^mu ) / denom;
                end
            end
        end
    end
    Y = U * ( M .* Atil ) * U';
end

function ZX = Zpinch(X)
    ZX = zeros(size(X));
    ZX(1:4,1:4) = X(1:4,1:4);
    ZX(5:8,5:8) = X(5:8,5:8);
end
