function global_stiff = assemble_matrix_const(dt, k_est_val, h , num_points,timescheme)
% k_est_val: estimated kernel values at dx:dx:ceil(delta/h)*h
    k = numel(k_est_val);
    k_est_val = [0; k_est_val];
    
    global_stiff = zeros(num_points, num_points);
    width = h;
    for i = 1:num_points
        for j = i - k:i + k
            if (j >= 1) && (j <= num_points) %&& (abs(j-i) >= 1)

                    global_stiff(i, j) = global_stiff(i, j) - k_est_val(abs(j-i)+1) * width; 
                    global_stiff(i, i) = global_stiff(i, i) + k_est_val(abs(j-i)+1) * width; 
            end
        end
    end
                    
    for i = 1:num_points
        if strcmp(timescheme,'BE')
            global_stiff(i,i) = global_stiff(i,i) + 1./dt/dt;
        elseif strcmp(timescheme, 'Newmark')
            global_stiff(i, i) = global_stiff(i,i) + 4. / dt / dt;
        end
    end

    % boundary condition
    for i  = 1:k
        global_stiff(i,:) = 0;
        global_stiff(i,i) = 1;
        global_stiff(num_points-i+1,:) = 0;
        global_stiff(num_points-i+1,num_points-i+1) = 1;
    end
end
