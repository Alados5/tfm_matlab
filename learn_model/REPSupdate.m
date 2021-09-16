function dw = REPSupdate(REWARDS)
% REWARDS is a row matrix with negative rewards (more negative -> worse)

options = optimset(@fmincon);
options.Display = 'off';
options.Algorithm = 'active-set';

% REPS
etamin = 0.0005;
etamax = 100;
etaini = 0.01;
Ekl = 0.5;

% Get weights with REPS
dualFunctionActual = @(eta_) dualfunction(eta_, REWARDS, Ekl);
eta2 = fmincon(dualFunctionActual, etaini, [], [], [], [], etamin, etamax,[],options);
if or(eta2==etamin, eta2==etamax)
    warning('Eta in its boundary')
end

dw = exp((REWARDS-max(REWARDS)*ones(size(REWARDS)))/eta2)';

end



