function loss = computeLoss(c,Abar,bbar,b_norm_sq)
% compute the loss for a given estimator c
 loss = c'*Abar*c - 2*c'*bbar + b_norm_sq; 
end
