function B = B_mat(inferInfo, SAVE_DIR)
% Output: B(i, j) = <phi_i, phi_j>_{L^2}
basis       = inferInfo.basis_funs; 
basis_knots = inferInfo.basis_funs_knots; 
basis_deg   = inferInfo.degree;
basis_str   = inferInfo.basis_str;
nBasis      = numel(basis);
delta       = inferInfo.supp_H(end);
filename    = [SAVE_DIR, 'B_mat_', basis_str, '.mat'];


if ~isfile(filename)
    B = zeros(nBasis);
            for i = 1:nBasis
                for j = i:nBasis
                    if contains(inferInfo.basis_str, 'Bspline')
                        xmin = min([basis_knots(i), basis_knots(j)]);
                        xmax = max([basis_knots(i+basis_deg), basis_knots(j+basis_deg)]);
                    elseif contains(inferInfo.basis_str, 'Bernstein')
                        xmin = 0;
                        xmax = delta;
                    end
                    fun = @(x) basis{i}(x) .* basis{j}(x);
                    %B(i,j) = inner_prod(basis{i}, basis{j}, [xmin xmax]);
                    B(i,j) = integral(fun, xmin, xmax);
                    B(j,i) = B(i,j);
                end
            end
    save(filename, 'B');
end
load(filename, 'B');

end

