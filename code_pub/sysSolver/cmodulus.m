function sum = cmodulus(mdelta,xi, k_est_val, dx)
% k_est_val(i) = K_est(xi)
c = k_est_val;

xx = abs(xi);
sum = 0.0d0;
for j=1:mdelta
    
  if(abs(xx-j*dx)<0.010d0*dx) 
    sum = sum+c(j);
  end
end

end