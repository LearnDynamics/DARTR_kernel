function [A_mat_new_rightPts, b_mat, B] = fast_update_matInfo(idx_nonzero_phi_conv, data, kernelInfo, obsInfo, inferInfo_old, inferInfo_new, SAVE_DIR)
    
    filename_matInfo_new = [SAVE_DIR, 'MatInfo_', kernelInfo.kernel_str, '_ui_', obsInfo.u_str,'_xj_', obsInfo.x_mesh_str, '_H_', inferInfo_new.basis_str, '.mat'];
    filename_matInfo_old = [SAVE_DIR, 'MatInfo_', kernelInfo.kernel_str, '_ui_', obsInfo.u_str,'_xj_', obsInfo.x_mesh_str, '_H_', inferInfo_old.basis_str, '.mat'];

    load(filename_matInfo_old); %, 'H_mat');    
    phi_old = inferInfo_old.basis_funs;
    phi_new = inferInfo_new.basis_funs;
    basis_knots = inferInfo_new.basis_funs_knots;
    basis_deg = inferInfo_new.degree;
    M_old = numel(phi_old); M_new = numel(phi_new);
    dist_x  = data.dist_x;
    dux_val = data.dux_val;
    dx = data.dx;
    A_mat_old_rightPts = A_mat_rightPts;
    A_mat_old_trapezoid = A_mat_trapezoid;

    B_old = B;
    [N, J, ~] = size(dux_val);
    
    A_mat_old_rightPts = reshape(A_mat_old_rightPts, J, N, M_old);
    A_mat_old_trapezoid = reshape(A_mat_old_trapezoid, J, N, M_old);
    A_mat_new_rightPts = zeros(J, N, M_new);
    A_mat_new_trapezoid = zeros(J, N, M_new);
    idx1 = idx_nonzero_phi_conv(1); idx2 = idx_nonzero_phi_conv(end);
    A_mat_new_rightPts(:, :, 1:idx1-1) = A_mat_old_rightPts(:, :, 1:idx1-1);
    A_mat_new_rightPts(:, :, M_new-M_old+idx2+1:end) = A_mat_old_rightPts(:, :, idx2+1:end);
    A_mat_new_trapezoid(:, :, 1:idx1-1) = A_mat_old_trapezoid(:, :, 1:idx1-1);
    A_mat_new_trapezoid(:, :, M_new-M_old+idx2+1:end) = A_mat_old_trapezoid(:, :, idx2+1:end);
    for n = 1:N
        parfor m = idx1:M_new-M_old+idx2
            a1 = dot(squeeze(phi_new{m}(dist_x)), squeeze(dux_val(n,:,:)), 1);
            A_mat_new_rightPts(:,n, m) = dx*a1;
            A_mat_new_trapezoid(:,n, m)= dx*(a1-0.5*a1(1)-0.5*a1(end));    
        end
    end  
    A_mat_new_rightPts = reshape(A_mat_new_rightPts, J*N, M_new);
    A_mat_new_trapezoid = reshape(A_mat_new_trapezoid, J*N, M_new);
    B = zeros(M_new, M_new);
    B(1:idx1-1, 1:idx1-1) = B_old(1:idx1-1, 1:idx1-1);
    B(M_new-M_old+idx2+1:end, M_new-M_old+idx2+1:end) = B_old(idx2+1:end, idx2+1:end);
    for i = 1:M_new
        for j = idx1:M_new-M_old+idx2
            xmin = min([basis_knots(i), basis_knots(j)]);
            xmax = max([basis_knots(i+basis_deg), basis_knots(j+basis_deg)]);
            B(i,j) = inner_prod(phi_new{i}, phi_new{j}, [xmin xmax]);
            B(j,i) = B(i,j);
        end
    end
    %{
    H_mat = zeros(J, M_new, N, J); 
    H_mat(:, 1:idx1-1, :) = H_mat_old(:, 1:idx1-1, :);
    H_mat(:, M_new-M_old+idx2+1:end, :) = H_mat_old(:, idx2+1:end, :);
    for n = 1:N
        parfor m = idx1:M_new-M_old+idx2
            H_mat(:,m,:) = squeeze(phi_new{m}(dist_x)).*squeeze(dux_val(n,:,:));
        end 
    end
    %}
    A_mat_rightPts  = A_mat_new_rightPts;
    A_mat_trapezoid = A_mat_new_trapezoid;
    save(filename_matInfo_new, 'A_mat_rightPts', 'A_mat_trapezoid', 'b_mat', 'B');%, 'H_mat');


end