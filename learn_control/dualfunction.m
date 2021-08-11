function [g] = dualfunction(eta, batch_return, epsilon)

n_batch = length(batch_return);
g = epsilon*eta+eta*(log(sum(exp((batch_return - max(batch_return))./eta))/n_batch)) +max(batch_return);

end