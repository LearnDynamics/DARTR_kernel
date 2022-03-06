function obsInfo = get_training_data(obsInfo, case_range, time_range)
%% get training data with case_range, time_range
ux_val  = obsInfo.u(case_range, time_range(2:end-1), :); 
utt_val = obsInfo.utt(case_range, time_range(2:end-1)-1, :); 
gx_val  = obsInfo.g(case_range, time_range(2:end-1), :); 
fx_val  = utt_val - gx_val;

%ux_val  = reshape(ux_val, N*(T-2), x_num);
%fx_val  = reshape(fx_val, N*(T-2), x_num);
%utt_val = reshape(utt_val, [], x_num);
data_str = obsInfo.training_data_str;
if numel(case_range) == 1
    d = 0;
else
    d = case_range(2)-case_range(1);
end
data_str = [data_str sprintf('case_%i_%i_%i_t_%i_%i_%i', case_range(1),d,case_range(end),time_range(1),time_range(2)-time_range(1),time_range(end))];
obsInfo.ux_val = ux_val;
obsInfo.fx_val = fx_val;
obsInfo.u_str = data_str;

obsInfo = rmfield(obsInfo, 'u');
obsInfo = rmfield(obsInfo, 'g');
obsInfo = rmfield(obsInfo, 'utt');
end
