function post_process_MOD_WET_EnKF_results(output_root_dir,...
    N_ensemble_members,output_state_maps,output_flux_maps,...
    output_flux_time_series,N_meas_periods,single_output_file_flag,...
    cleanup_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to post-process the raw EnKF outputs files. The raw 
% files were stored individually for each ensemble member over each 
% measurement period during application of the EnKF. This code merges the 
% outputs across measurement periods and then either stores either 
% replicate-wise or as one large output file.

% Inputs:
% output_root_dir: directory where raw outputs are stored
% N_ensemble_members: number of ensemble members
% output_state_maps: output states to merge
% output_flux_maps: output fluxes to merge
% output_flux_time_series: output flux time series to merge
% N_meas_periods: number of measurement periods
% single_output_file_flag: Set equal to 0 to store all merged outputs in 
%    one file (will be a large file). Set equal to 1 to store merged 
%   outputs in individual files for each replicate.
% cleanup_flag: Set equal to 0 to leave raw data files or 1 to delete them

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load and merge outputs
% Loop over ensemble
for i_ensemble=1:N_ensemble_members

    disp(['Post-processing ensemble ' num2str(i_ensemble) ' out of ' ...
        num2str(N_ensemble_members)])
    % Loop over measurement periods
    for i_meas_period=1:N_meas_periods
            % Name of stored output file for this measurement period/ensemble member
            filename=[output_root_dir 'ensemble_member_' ...
                num2str(i_ensemble) '_meas_period_' ...
                num2str(i_meas_period) '_output.mat'];
            load(filename)
            
            if (i_meas_period==1)
                ind1=1;
            else
                ind1=ind2+1;
            end
            
            command=['ind2=ind1+size(states.maps.' ...
                output_state_maps(1).names ',3) - 1;'];
            eval(command);
            
            if (single_output_file_flag==0)
                index=1;
            else
                index=i_ensemble;
            end
            % Grab state maps and store them
            for istate=1:length(output_state_maps)
                command=['STATE_maps.' output_state_maps(istate).names ...
                    '(:,:,ind1:ind2,index)=states.maps.' ...
                    output_state_maps(istate).names ';'];
                eval(command)
            end      
            % Grab flux maps and store them
            for iflux=1:length(output_flux_maps)
                command=['FLUX_maps.' output_flux_maps(iflux).names ...
                    '(:,:,ind1:ind2,index)=fluxes.maps.' ...
                    output_flux_maps(iflux).names ';'];
                eval(command)
            end
            % Grab hydrograph time series and store it
            for iflux=1:length(output_flux_time_series)
                command=['FLUX_time_series.' ...
                    output_flux_time_series(iflux).names ...
                    '(ind1:ind2,index)=fluxes.time_series.' ...
                    output_flux_time_series(iflux).names ';'];
                eval(command)
            end
            
    end % end measurement interval loop
    
    % Store individual replicate files (depending on flag)
    if (single_output_file_flag==0)
        % Store merged output file for each replicate
        output_filename=[output_root_dir 'ensemble_member_' ...
            num2str(i_ensemble) '_merged_output.mat'];
        save(output_filename,'STATE_maps','FLUX_maps','FLUX_time_series')
    end
end % end ensemble replicate loop

% Store single file (depending on flag)
if (single_output_file_flag==1)
    output_filename=[output_root_dir 'merged_EnKF_outputs.mat'];
    save(output_filename,'STATE_maps','FLUX_maps','FLUX_time_series')
end

% Cleanup raw files (depending on flag)
if (cleanup_flag==1)
    for i_ensemble=1:N_ensemble_members
        % Loop over measurement periods
        for i_meas_period=1:N_meas_periods
            % Name of stored output file for this measurement period/ensemble member
            filename=[output_root_dir 'ensemble_member_' ...
                num2str(i_ensemble) '_meas_period_' ...
                num2str(i_meas_period) '_output.mat'];
            delete(filename);
        end
    end
end

return