function [u_normed, f_normed] = normalize_byL2norm(u, f, dx)
d = numel(size(u));
n = vecnorm(u, 2, d)*sqrt(dx);
u_normed = u./n;  
f_normed = f./n;

end