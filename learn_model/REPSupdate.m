function dw=REPSupdate(REWARDS)

% REWARDS is a row matrix with negative rewards (more negative -> worse)
Ekl=0.5;

options = optimset(@fmincon);
options.Display='off';
options.Algorithm='active-set';

%% REPS

etamin = 0.0005;
etamax = 100;
etaini = 0.01;

% get weights with REPS
dualFunctionActual = @(eta_) dualfunction(eta_, REWARDS, Ekl);
eta2 = fmincon(dualFunctionActual, etaini, [], [], [], [], etamin, etamax,[],options);
if or(eta2==etamin,eta2==etamax)
    warning('Eta in its boundary')
    % eta2
end
dw = exp((REWARDS-max(REWARDS)*ones(size(REWARDS)))/eta2)';
% dw -> relative weights

% Z = (sum(dw)*sum(dw) - sum(dw .^ 2))/sum(dw);
% mw = sum(bsxfun(@times, WEIGHTS', dw)',2)./sum(dw);
% summ=0;
% for ak=1:size(WEIGHTS,2)
%     summ=summ+dw(ak)*((WEIGHTS(:,ak)-mw)*(WEIGHTS(:,ak)-mw)');%/sum(dw);
% end
% Sw = summ./(Z+1e-9);

end



