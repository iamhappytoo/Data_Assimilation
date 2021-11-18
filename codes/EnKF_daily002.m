function MOD_WET_wrapper_with_EnKF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Steve Margulis, May 5, 2017
%
% This is a "wrapper" script that allows for application of the EnKF with
% the MOD-WET watershed model. The code is setup to be adapted to include 
% various uncertainties in initial conditions, time-invariant parameters, 
% and meteorological forcings. 
%
% For tractability, the current setup uses lognormal multiplicative 
% perturbations as a way of adding uncertainty to model inputs. 
% Modifications could be made to include additive errors of varying form.
% The initial condition and parameters uncertainty is introduced in the 
% form of spatially-correlated fields (that can be cross-correlated) using
% the TBCOSIM function (tbcosim_uncond_func.m) that has been adapted for 
% simulation of random 
%
% Calling of MOD-WET (including re-starts with updated states) is handled
% entirely through modification of its input structure array: restart_info.
% It is assumed that MOD-WET is already setup for a nominal watershed with
% a nominal set of parameters and forcings. Those nominal inputs are first 
% loaded, then perturbations are generated and those inputs that are 
% changed are saved to the restart_info array for initial conditions and 
% parameters and to meteorological forcing input files.
% 
% The EnKF code is embedded at the end of the main function and it is 
% assumed a measurement_model.m function is defined with the appropriate 
% measurement model. The measurements (along with the appropriate time
% indices, measurement error, etc.) are loaded. For synthetic experiments,
% this file can be generated from the MOD_WET_wrapper_truegen.m function.
%
% For simplicity all raw outputs and perturbed inputs are saved to files 
% for later use. A post-processinsg code is called at the end to merge the
% results across measurement (i.e. propagation) intervals for each
% replicate. A flag allows for storing the desired states/fluxes in
% individual files for each replicate or one (potentially large) file with
% all replicates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;close all
%% Specify path to MOD-WET functions
% Specify  to root directory for the MOD-WET functions
toolbox_path='C:/Users/BANZH/Downloads/2019spring/CEE steve/MOD_WET_2017a/MOD_WET_2017a/MOD_WET_src';
% Add path
addpath(genpath(toolbox_path))
% Specify root directory where to store outputs
output_root_dir='C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/long_enkf_output_daily_0.02/';
mkdir(output_root_dir);
%% --------- SPECIFY INPUTS HERE ----------------------------------------%%
% Note: Additional inputs are specified in: model_input_uncertainty_parameters.m 
% and in the measurement model: measurement_model.m
%
% Specify number of ensemble members to use with EnKF
N_ensemble_members=20;

% Localization flag (set to 1 to perform localization; 0 to skip
% localization). Note: The localization used herein assumes an exponential
% spatial correlation function.
localization_parameters.localization_flag=0;
if (localization_parameters.localization_flag==1)
    % Set localization range parameter. This is the exponential correlation
    % length parameter in pixel distance. The number of pixels is converted 
    % to an actual distance based on the pixel dimensions below during the
    % localization application.
    localization_range_parameter=5; 
end

%% Load measurements and measurement info. (i.e. from true simulation)
measurement_input_file='C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/long_true_output_daily_0.02/measurements.mat';
% This loads: meas_day_vector, ZMEAS, Cv
load(measurement_input_file, ...
    'meas_day_vector', 'ZMEAS','Cv','start_day','n_days')

day_start_vector=[0 meas_day_vector n_days]+start_day;
N_meas_periods=length(day_start_vector)-1;

%% Define all MOD-WET states
mod_wet_states(1).names='Srz';mod_wet_states(2).names='Suz';
mod_wet_states(3).names='SD';mod_wet_states(4).names='Tsurf';
mod_wet_states(5).names='SWE';mod_wet_states(6).names='snowdepth';
mod_wet_states(7).names='snowfrac';mod_wet_states(8).names='snowdens';
mod_wet_states(9).names='Td';
mod_wet_states(10).names='NDayLastSnow';

%% Define uncertain model inputs. Uncertainty parameters are specified in the function: model_input_uncertainty_parameters.m
% This is the subset of states that are treated as uncertain and that will
% be updated via the DA, e.g.:
% Updating of Srz and SD only:
uncertain_state_numbers=[1,3]; % here Srz and SD from list above.
% Updating SWE only:
%uncertain_state_numbers=[5]; % here SWE from list above.


for ii=1:length(uncertain_state_numbers)
    uncertain_states(ii).names=mod_wet_states((uncertain_state_numbers(ii))).names;
end
nvar=length(uncertain_states);
%
% Define uncertain time-invariant parameters:
% These are not updated, but add uncertainty to the states via the model
% physics
% these two parameters are really coupled, but can control this by making them perfectly correlated below
uncertain_params(1).names='T0';
uncertain_params(2).names='K0'; 
%
% Define uncertain time-varying forcing:
% Meteorological inputs that are treated as uncertain. Note: As currently
% constructed MOD-WET takes in only a single set of forcings and then
% distributed the forcing spatially within the model. This makes it
% difficult to treat these inputs as random spatial fields. SO instead they
% are just treated as correlated single perturbations that can vary in
% time, but are deterministically distributed in space according to the 
% MOD-WET disaggregation equations.
uncertain_forcings(1).names='PPT';

%% Specify states/fluxes to post-process and save from the EnKF results.
% Specify state maps to post-process: 
% Note: These are not necessarily only the states that were included in the 
% EnKF state vector, but could include those that were not explicitly 
% updated, but were implicitly updated via the model physics.
output_state_maps(1).names='Srz';output_state_maps(2).names='SD'; output_state_maps(3).names='SWE';

% Specify flux maps:
% Note: These variables are all implicitly updated via the state updates.
% Often these are of equal, if not more, interest than states.
output_flux_maps(1).names='ET';

% Specify time series (here only outlet hydrograph)
output_flux_time_series(1).names='outlet_hydrograph';

% Specify flag whether to store all merged outputs as one file per 
% measurement interval or as one large file.
single_output_file_flag=1; % 0=ensemble-wise files; 1= one single file. 

% Specify flag whether to cleanup the raw outputs and keep only the merged
% output files.
cleanup_flag=0; % 0=keep all files; 1=delete raw output files
%% --------- END OF INPUT SPECIFICATION ---------------------------------%%

%% Call initialization codes to get nominal parameters for model simulation. 
% This first part is just a simple mechanism to load the nominal 
% parameters/initial from the existing MOD-WET parameter/initial condition 
% functions. The desired variable values should be set in those functions 
% and are loaded here. By running these functions all of the parameters 
% into the "restart_info" structure array. The nominal values can then be 
% perturbed and replaced as inputs. 
%
% Load nominal control parameters (i.e. start day, time step, number of 
% days, etc.)
[restart_info.control_params]=set_control_parameters;
% This is to make sure hourly maps are output (assumes forcing is provided
% at 0.25 hour interval)
restart_info.control_params.map_frq2store=4;
%
% Load nominal static physical parameters (soil parameters, TOPMODEL 
% parameters, etc. that are defined in the standard MOD-WET input
% functions.
[restart_info.params]=set_static_physical_parameters(...
    restart_info.control_params);
%
% pre-allocation of variables
[restart_info.states,~,~,~]= ...
                pre_allocate_variables(restart_info.control_params,...
                restart_info.params);
%
% Load nominal initial conditions for state variables (i.e. "old_vars")
[~,~,restart_info.old_vars]= ...
        initialize_model(restart_info.control_params,restart_info.params,...
        restart_info.states);

%% Grab domain-specific parameters
% Grab size of domain for spatial field perturbations
%[ny,nx]=size(restart_info.params.static_maps.mask);
% Grab mask of pixel indices (the length of this vector is the size of 
% the state vector per variable.
pixel_indices=find(restart_info.params.static_maps.mask==1);
Npix=length(pixel_indices);

%% Pre-allocate state/measurement vectors
% state vector
Y_ENS=NaN(Npix*length(uncertain_states),N_ensemble_members);
% meas. vector (this will depend on what the measurements is
Z_ENS=NaN(5,N_ensemble_members);

%% Compute some inputs needed for localization
if (localization_parameters.localization_flag==1)
    % Coordinate arrays for domain
    [EASTING,NORTHING]=meshgrid(restart_info.params.static_maps.easting,...
        restart_info.params.static_maps.northing);
    dx=restart_info.params.static_maps.easting(2)-...
        restart_info.params.static_maps.easting(1);
    % Define localization correlation function (i.e. using exponential
    % correlation function)
    % in meters (dx=resolution of DEM)
    range_parameter=localization_range_parameter*dx; 
    
    % Coordinates (easting/northing) for model pixels
    X=[EASTING(pixel_indices) NORTHING(pixel_indices)]';

    % Note: Since the state vector includes two states at each pixel we 
    % need to include pixel coordinates twice; this is simply the 
    % coordinates of every element of the state vector.
    XX=repmat(X,[1 length(uncertain_state_numbers)]);
    
    % calculate a distance matrix between all coordinate pairs in coordinate
    % vector. This is just an efficient distance calculator using dot products.
    %distance_yy = sqrt(bsxfun(@plus,dot(XX,XX,1)',dot(XX,XX,1))-2*(XX'*XX));
    %localization_parameters.rho_yy_localization=exp(-distance_yy/range_parameter);
    
    % Coordinates of measurements
    ZZ=[X];
    % calculate distance matrix between state locations and measurement
    % locations
    distance_yz = sqrt(bsxfun(@plus,dot(XX,XX,1)',dot(ZZ,ZZ,1))-2*(XX'*ZZ));
    localization_parameters.rho_yz=exp(-distance_yz/range_parameter);
    
    % calculate distance matrix between measurement locations
    distance_zz = sqrt(bsxfun(@plus,dot(ZZ,ZZ,1)',dot(ZZ,ZZ,1))-2*(ZZ'*ZZ));
    localization_parameters.rho_zz=exp(-distance_zz/range_parameter);

    clear X XX distance_yy range_parameter ZZ distance_yz distance_zz
end

%% MAIN LOOPING -- ACROSS MEAS. INTERVALS AND FULL ENSEMBLE
%% Loop through measurement intervals
for i_meas_period=1:N_meas_periods
      
    %% Loop through an ensemble of simulations (with different inputs)        
    for i_ensemble=1:N_ensemble_members
         %% Call initialization codes to get nominal parameters for model simulation. 
        % This first part is just a simple mechanism to load the nominal 
        % parameters/initial from the existing MOD-WET parameter/initial condition 
        % functions. The desired variable values should be set in those functions 
        % and are loaded here. By running these functions all of the parameters 
        % into the "restart_info" structure array. The nominal values can then be 
        % perturbed and replaced as inputs. 
        %
        % Load nominal control parameters (i.e. start day, time step, number of 
        % days, etc.)
        [restart_info.control_params]=set_control_parameters;
        % This is to make sure hourly maps are output (assumes forcing is provided
        % at 0.25 hour interval)
        restart_info.control_params.map_frq2store=4;
        %
        % Load nominal static physical parameters (soil parameters, TOPMODEL 
        % parameters, etc. that are defined in the standard MOD-WET input
        % functions.
        [restart_info.params]=set_static_physical_parameters(...
            restart_info.control_params);
        %
        % pre-allocation of variables
        [restart_info.states,~,~,~]= ...
                        pre_allocate_variables(restart_info.control_params,...
                        restart_info.params);
        %
        % Load nominal initial conditions for state variables (i.e. "old_vars")
        [~,~,restart_info.old_vars]= ...
                initialize_model(restart_info.control_params,restart_info.params,...
                restart_info.states);

        % Set period over which to run MOD-WET for the propagation step
        start_day=day_start_vector(i_meas_period);
        n_days=day_start_vector(i_meas_period+1)-start_day;

        %% Set simulation control parameters:
        % Replace the starting day of the simulation (0 is the beg. of the WY)
        restart_info.control_params.start_day=start_day; % start simulation on Nov. 1st 
        % Replace the length of the simulation and related input parameters.
        restart_info.control_params.n_days=n_days; % through end of Dec.
        restart_info.control_params.ntime=1./restart_info.control_params.dt*24.*restart_info.control_params.n_days;
        % Length of time series to store (needs to be an integer)
        restart_info.control_params.nt=restart_info.control_params.ntime/restart_info.control_params.timeseries_frq2store; 
        % Set index to initialize forcings 
        restart_info.control_params.start_time=1+...     
        restart_info.control_params.start_day*24./restart_info.control_params.dt;

    % Define nominal forcing filename (for re-using below)
    nominal_forcing_filename=restart_info.control_params.met_data_filename;
        % Prior to start of assimilation, perturb. uncertain inputs (only needs to be done once):
        if (i_meas_period==1) 

            % Replace output filename with one specific to this ensemble member
            restart_info.control_params.output_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_meas_period_' num2str(i_meas_period) '_output.mat'];

            %% Perturb the  initial conditions, physical parameters, and forcing:            
            %
            % load nominal forcings
            eval(['load ' '''' nominal_forcing_filename ''''])

            %% Call function that generates perturbations based on specified statistics
            seed1=i_ensemble;
            seed2=i_ensemble+N_ensemble_members;
            [IC_perturbations,params_perturbations,forcing_perturbations]=...    
            model_input_uncertainty_parameters(uncertain_states,uncertain_params,...
            uncertain_forcings,restart_info,length(PPT),seed1,seed2);
            %
            % ---- Perturb initial conditions (based on statistics above):
            % Loop through states that will be perturbed
            for ii=1:length(uncertain_states)
                command=['restart_info.old_vars.' uncertain_states(ii).names '0 = IC_perturbations(:,:,ii).* restart_info.old_vars.' uncertain_states(ii).names '0;'];
                eval(command)
            end
            %% NEW CODE !!
            % Make sure any updated states are within defined (physical)
            % bounds
            % Srz (needs to be within min-max range)
            I=find(restart_info.old_vars.Srz0<restart_info.params.static_maps.Srzmin);
            restart_info.old_vars.Srz0(I)=restart_info.params.static_maps.Srzmin(I);
            I=find(restart_info.old_vars.Srz0>restart_info.params.static_maps.Srzmax);
            restart_info.old_vars.Srz0(I)=restart_info.params.static_maps.Srzmax(I);
            % SD (should be non-negative)
            I=find(restart_info.old_vars.SD0<0);
            restart_info.old_vars.SD0(I)=0;
            % Set Suz0 to zero where SD0 is zero
            restart_info.old_vars.Suz0(I)=0;
            % SWE and other snow states (should be non-negative)
            I=find(restart_info.old_vars.SWE0<0);
            restart_info.old_vars.SWE0(I)=0;
            restart_info.old_vars.snowdens0(I)=0;
            restart_info.old_vars.snowdepth0(I)=0;
            restart_info.old_vars.snowfrac0(I)=0;
            % Make sure temperature in snow-covered areas is less than or
            % equal to freezing
            I=find(restart_info.old_vars.SWE0>0 & restart_info.old_vars.Tsurf0>273.15);
            restart_info.old_vars.Tsurf0(I)=273.15;
            %% END NEW CODE !!
            %
            % ---- Perturb time-invariant parameters (based on statistics above):
            % Loop through parameters that will be perturbed
            for ii=1:length(uncertain_params)
                command=['restart_info.params.static_maps.' uncertain_params(ii).names ' = params_perturbations(:,:,ii).* restart_info.params.static_maps.' uncertain_params(ii).names ';'];
                eval(command)
            end
            %% NEW CODE !! 
            % save parameters (within restart_info structure array) to file
            % for later use
            restart_params_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_perturbed_params.mat'];
            perturbed_params=restart_info.params;
            eval(['save ' '''' restart_params_filename '''' ' perturbed_params'])
            %% END NEW CODE !!
            
            %
            % ---- Perturb time-varying forcings (based on statistics above):
            % Loop through forcings and apply perturbations
            for ii=1:length(uncertain_forcings)
                command=[uncertain_forcings(ii).names ' = ' uncertain_forcings(ii).names '.* forcing_perturbations(ii,:);' ];
                eval(command)
            end
            %
            % Write forcings to file that will be loaded by MOD-WET
            restart_info.control_params.met_data_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_met_forcing_input.mat'];
            eval(['save ' '''' restart_info.control_params.met_data_filename '''' ' time gage_elev PPT SW Ta qa U Psfc'])

        else % otherwise use already perturbed values or updated states at later propagation steps

            %% NEW CODE !! 
            % Load perturbed parameters (within restart_info structure array)
            restart_params_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_perturbed_params.mat'];
            eval(['load ' '''' restart_params_filename '''' ' perturbed_params'])
            restart_info.params=perturbed_params;
            %% END NEW CODE !!            
            
            % Set proper forcing input file for this ensemble member
            restart_info.control_params.met_data_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_met_forcing_input.mat'];
            
            % Load old output from previous measurement period to
            % re-initialize non-updated states
            restart_info.control_params.output_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_meas_period_' num2str(i_meas_period-1) '_output.mat'];
            % Load outputs from this simulation to initialize states
            eval(['load ' '''' restart_info.control_params.output_filename '''' ' states'])
            %
            % Re-initialize states with values at end of previous period
            for ii=1:length(mod_wet_states)    
                if (ii==5 | ii==6) % SWE or snow depth, convert to mm for input to model:
                    command=['restart_info.old_vars.' mod_wet_states(ii).names '0 = states.maps.' mod_wet_states(ii).names '(:,:,end)*1000.;'];
                elseif (ii==10) % For the number of days since albedo was re-set; this has no trailing "0" in the variable in MOD-WET
                    command=['restart_info.old_vars.' mod_wet_states(ii).names ' = states.maps.' mod_wet_states(ii).names '(:,:,end);'];
                else
                    command=['restart_info.old_vars.' mod_wet_states(ii).names '0 = states.maps.' mod_wet_states(ii).names '(:,:,end);'];
                end
                eval(command)
            end
            
            % Overwrite old ICs with updated states
            Y_IC=Y_UPDATE(:,i_ensemble);
            Y_IC=reshape(Y_IC,Npix,nvar);            
            for ii=1:length(uncertain_states)
                tmp=NaN(size(restart_info.params.static_maps.mask));
                tmp(pixel_indices)=Y_IC(:,ii);
                command=['restart_info.old_vars.' uncertain_states(ii).names '0 = tmp;'];
                eval(command)
            end

            %% NEW CODE !!
            % Make sure any updated states are within defined (physical)
            % bounds
            % Srz (needs to be within min-max range)
            I=find(restart_info.old_vars.Srz0<restart_info.params.static_maps.Srzmin);
            restart_info.old_vars.Srz0(I)=restart_info.params.static_maps.Srzmin(I);
            I=find(restart_info.old_vars.Srz0>restart_info.params.static_maps.Srzmax);
            restart_info.old_vars.Srz0(I)=restart_info.params.static_maps.Srzmax(I);
            % SD (should be non-negative)
            I=find(restart_info.old_vars.SD0<0);
            restart_info.old_vars.SD0(I)=0;
            % Set Suz0 to zero where SD0 is zero
            restart_info.old_vars.Suz0(I)=0;
            % SWE and other snow states (should be non-negative)
            I=find(restart_info.old_vars.SWE0<0);
            restart_info.old_vars.SWE0(I)=0;
            restart_info.old_vars.snowdens0(I)=0;
            restart_info.old_vars.snowdepth0(I)=0;
            restart_info.old_vars.snowfrac0(I)=0;
            % Make sure temperature in snow-covered areas is less than or
            % equal to freezing
            I=find(restart_info.old_vars.SWE0>0 & restart_info.old_vars.Tsurf0>273.15);
            restart_info.old_vars.Tsurf0(I)=273.15;
            %% END NEW CODE !!
            
            % Create output filename with one specific to this ensemble member
            restart_info.control_params.output_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_meas_period_' num2str(i_meas_period) '_output.mat'];

        end
        
        %% Call MOD-WET with perturbed inputs
        MOD_WET_model_driver(restart_info)

        %% Load model outputs
        eval(['load ' '''' restart_info.control_params.output_filename '''' ' states'])

        %% Collect state vector for update time
        % Create state vector
        Y=NaN(Npix,length(uncertain_states));
        for ii=1:length(uncertain_states)
            % grab states at end of propagation step to be updated
            if (strcmp(uncertain_states(ii).names,'SWE') || strcmp(uncertain_states(ii).names,'snowdepth'))
                command=['tmp=states.maps.' uncertain_states(ii).names '(:,:,end)*1000;']; % converts from meters to mm
            else
                command=['tmp=states.maps.' uncertain_states(ii).names '(:,:,end);'];
            end
            eval(command)
            command=['Y(:,ii)=tmp(pixel_indices);'];
            eval(command)
        end
        Y_ENS(:,i_ensemble)=Y(:); % store as a single column vector
        
        %% Collect predicted measurement vector for update time
        Z_ENS(:,i_ensemble)=measurement_model(Y);  
        
    end
    %% end ensemble loop

    %% Call EnKF (only done when measurements exist at end of period)
    if (i_meas_period<N_meas_periods)
        % Grab measurement at this time
        zmeas=ZMEAS(:,i_meas_period);
        % Update states via EnKF
        [Y_UPDATE]=ENKF(Y_ENS,Z_ENS,zmeas,Cv,localization_parameters);
        %pause
        %'Y_ENS'
        %Y_ENS
        %'Z_ENS'
        %Z_ENS
        %pause
        %'Y_UPDATE'
        %Y_UPDATE
        %pause
    end
end
%% end measurement period loop

%% Call post-processing function to merge outputs
post_process_MOD_WET_EnKF_results(output_root_dir,...
    N_ensemble_members,output_state_maps,output_flux_maps,...
    output_flux_time_series,N_meas_periods,single_output_file_flag,...
    cleanup_flag)

return

%% Embedded functions
function [Y_UPDATE]=ENKF(Y_ENS,Z_ENS,ZMEAS,Cv,localization_parameters)
    % Apply the EnKF to a generic ensemble of states and predicted measurement
    
    % Grab size of ensemble
    Nreps=size(Y_ENS,2);
    % Pre-allocate update 
    Y_UPDATE=NaN(size(Y_ENS));
    % Define measurement error covariance assuming homoscedastic error with 
    % provided Cv as the diagonal (variance)
    %Cv=(0.02)^2
    CV=diag(Cv*ones(size(Z_ENS,1),1));

    % Compute sample mean vectors
    ymean=mean(Y_ENS,2)*ones(1,Nreps);
    zmean=mean(Z_ENS,2)*ones(1,Nreps);
    % Compute sample covariance matrices
    Cyz=((Y_ENS-ymean)*(Z_ENS-zmean)')/(Nreps-1);
    Czz=((Z_ENS-zmean)*(Z_ENS-zmean)')/(Nreps-1);
    
    % If performing localization:
    if (localization_parameters.localization_flag==1)
        Cyz=localization_parameters.rho_yz.*Cyz;
        Czz=localization_parameters.rho_zz.*Czz;
    end
    
    % compute Kalman gain
    kalm=Cyz/(Czz+CV);
    
    % Measurement error realizations needed for update
    v=mvnrnd(zeros(1,size(Z_ENS,1)),CV,Nreps)';
    
    % Update states
    for k=1:Nreps
        Y_UPDATE(:,k)=Y_ENS(:,k)+kalm*(ZMEAS+v(:,k)-Z_ENS(:,k));
    end

return


%% Plot and compare the updated results with truth
load('C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/long_enkf_output_daily_0.02/merged_EnKF_outputs.mat');
EnsSD=nan(26,19,5832,20);
EnsSrz=nan(26,19,5832,20);
EnsET=nan(26,19,5832,20);
EnsQ=nan(5832,20);
EnsSD=STATE_maps.SD;
EnsSrz=STATE_maps.Srz;
EnsET=FLUX_maps.ET;
EnsQ=FLUX_time_series.outlet_hydrograph;
MeanEnsSD=mean(EnsSD,4);
MeanEnsSrz=mean(EnsSrz,4);
MeanEnsET=mean(EnsET,4);
MeanEnsQ=mean(EnsQ,2);
load('C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/long_true_output_daily_0.02/measurements.mat');
load('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/long_true_output_daily_0.02/true_output.mat');
Srz2=(EnsSrz).^2;
SD2=(EnsSD).^2;
ET2=(EnsET).^2;
Q2=(EnsQ).^2;
Srz_stddev=sqrt(abs(mean(Srz2,4)-MeanEnsSrz.^2));
SD_stddev=sqrt(mean(SD2,4)-MeanEnsSD.^2);
ET_stddev=sqrt(abs(mean(ET2,4)-MeanEnsET.^2));
Q_stddev=sqrt(abs(mean(Q2,2)-MeanEnsQ.^2));
timeind=42;
figure(1);
subplot(1,2,1);
imagesc(states.maps.SD(:,:,timeind*24));axis image;
colorbar;title(['True SD at end of Day' num2str(timeind)])
subplot(1,2,2);imagesc(MeanEnsSD(:,:,timeind*24));axis image;
colorbar;title(['Ensemble mean SD at end of Day' num2str(timeind)])
pause
figure(2);
subplot(1,2,1);
imagesc(states.maps.Srz(:,:,timeind*24));axis image;
colorbar;title(['True Srz at end of Day' num2str(timeind)])
subplot(1,2,2);imagesc(MeanEnsSrz(:,:,timeind*24));axis image; 
colorbar;title(['Ensemble mean Srz at end of Day' num2str(timeind)])
pause
figure(3);
subplot(1,2,1);
imagesc(fluxes.maps.ET(:,:,timeind*24));axis image;
colorbar;title(['True ET at end of Day' num2str(timeind)])
subplot(1,2,2);imagesc(MeanEnsET(:,:,timeind*24));axis image; 
colorbar;title(['Ensemble mean ET at end of Day' num2str(timeind)])
pause 

DOY=1:1:5832
row0=19;
col0=3;
% plot at a pixel that have SD measurements, for [SD]
figure(4)
plot(DOY,squeeze(states.maps.SD(row0,col0,:)),'b','LineWidth',col0);hold on
plot(meas_day_vector*24,ZMEAS(1,:),'mo','LineWidth',2,'MarkerSize',10);hold on
plot(DOY,squeeze(MeanEnsSD(row0,col0,:)),'k','LineWidth',col0);hold on
plot(squeeze(MeanEnsSD(row0,col0,:)+2*SD_stddev(row0,col0,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);
plot(squeeze(MeanEnsSD(row0,col0,:)-2*SD_stddev(row0,col0,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1); 
xlabel('Model time step');ylabel('SD (m)')
title(['Time series and measurements for pixel: (' num2str(row0) ',' num2str(col0) ')'])
set(gca,'FontSize',14);grid
legend('True','Meas.','EnKF Ens. Mean','+/- 2 Std. Dev.')

% plot at a pixel that have SD measurements, for [Srz]
figure(5)
plot(DOY,squeeze(states.maps.Srz(row0,col0,:)),'b','LineWidth',3);hold on
plot(DOY,squeeze(MeanEnsSrz(row0,col0,:)),'k','LineWidth',3);hold on
plot(squeeze(MeanEnsSrz(row0,col0,:)+2*Srz_stddev(row0,col0,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);
plot(squeeze(MeanEnsSrz(row0,col0,:)-2*Srz_stddev(row0,col0,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1); 
xlabel('Model time step');ylabel('Srz (m)')
title(['Time series and measurements for pixel: (' num2str(row0) ',' num2str(col0) ')'])
set(gca,'FontSize',14);grid
legend('True','EnKF Ens. Mean','+/- 2 Std. Dev.')

% plot at a pixel that have SD measurement, for [ET]
figure(6)
plot(DOY,squeeze(fluxes.maps.ET(row0,col0,:)),'b','LineWidth',1);hold on
plot(DOY,squeeze(MeanEnsET(row0,col0,:)),'k','LineWidth',1);hold on
plot(squeeze(MeanEnsET(row0,col0,:)+2*ET_stddev(row0,col0,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);
plot(squeeze(MeanEnsET(row0,col0,:)-2*ET_stddev(row0,col0,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1); 
xlabel('Model time step');ylabel('ET (m/h)')
title(['Time series and measurements for pixel: (' num2str(row0) ',' num2str(col0) ')'])
set(gca,'FontSize',14);grid
legend('True','EnKF Ens. Mean','+/- 2 Std. Dev.')
row=15;
col=15;
% plot some pixels that are not measured [SD]
figure(7)
plot(DOY,squeeze(states.maps.SD(row,col,:)),'b','LineWidth',3);hold on
%plot(meas_day_vector*24,ZMEAS(1,:),'mo','LineWidth',2,'MarkerSize',10);hold on
plot(DOY,squeeze(MeanEnsSD(row,col,:)),'k','LineWidth',3);hold on
plot(squeeze(MeanEnsSD(row,col,:)+2*SD_stddev(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);hold on
plot(squeeze(MeanEnsSD(row,col,:)-2*SD_stddev(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);hold on
for i=1:6
plot(squeeze(EnsSD(row,col,:,i)),'c');hold on
end
xlabel('Model time step');ylabel('SD (m)')
title(['Time series and measurements for pixel: (' num2str(row) ',' num2str(col) ')'])
set(gca,'FontSize',14);grid
legend('True','EnKF Ens. Mean','+2 Std. Dev.','-2 Std. Dev.','Ens')
% plot some pixels that are not measured [Srz]
figure(8)
plot(DOY,squeeze(states.maps.Srz(row,col,:)),'b','LineWidth',3);hold on
%plot(meas_day_vector*24,ZMEAS(1,:),'mo','LineWidth',2,'MarkerSize',10);hold on
plot(DOY,squeeze(MeanEnsSrz(row,col,:)),'k','LineWidth',3);hold on
plot(squeeze(MeanEnsSrz(row,col,:)+2*Srz_stddev(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);hold on
plot(squeeze(MeanEnsSrz(row,col,:)-2*Srz_stddev(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);hold on
for i=1:6
plot(squeeze(EnsSrz(row,col,:,i)),'c');hold on
end
xlabel('Model time step');ylabel('Srz (m)')
title(['Time series and measurements for pixel: (' num2str(row) ',' num2str(col) ')'])
set(gca,'FontSize',14);grid
legend('True','EnKF Ens. Mean','+2 Std. Dev.','-2 Std. Dev.','Ens')

% plot some pixels that are not measured [ET]
figure(9)
plot(DOY,squeeze(fluxes.maps.ET(row,col,:)),'b','LineWidth',3);hold on
%plot(meas_day_vector*24,ZMEAS(1,:),'mo','LineWidth',2,'MarkerSize',10);hold on
plot(DOY,squeeze(MeanEnsET(row,col,:)),'k','LineWidth',3);hold on
plot(squeeze(MeanEnsET(row,col,:)+2*ET_stddev(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);hold on
plot(squeeze(MeanEnsET(row,col,:)-2*ET_stddev(row,col,:)),'--','Color',[0.85 0.33 0.1],'LineWidth',1);hold on
for i=1:6
plot(squeeze(EnsET(row,col,:,i)),'c');hold on
end
xlabel('Model time step');ylabel('ET (m/h)')
title(['Time series and measurements for pixel: (' num2str(row) ',' num2str(col) ')'])
set(gca,'FontSize',14);grid
legend('True','EnKF Ens. Mean','+2 Std. Dev.','-2 Std. Dev.','Ens')

figure(10)
plot(DOY,squeeze(fluxes.time_series.outlet_hydrograph),'b','LineWidth',3);hold on
plot(DOY,squeeze(MeanEnsQ),'k','LineWidth',3);hold on
plot(squeeze(MeanEnsQ+3*Q_stddev),'--','Color',[0.85 0.33 0.1],'LineWidth',1);hold on
plot(squeeze(MeanEnsQ-3*Q_stddev),'--','Color',[0.85 0.33 0.1],'LineWidth',1);hold on
for i=1:6
plot(squeeze(EnsQ(:,i)),'c');hold on
end
xlabel('Model time step');ylabel('Q (m^3/s)')
title(['Time series of Q at outlet'])
set(gca,'FontSize',14);grid
legend('True','EnKF Ens. Mean','+3 Std. Dev.','-3 Std. Dev.','Ens')

cd('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\long_enkf_output_daily_0.02\');
mkdir('C:\Users\BANZH\Downloads\2019spring\CEE Steve\Project_outputs\long_enkf_output_daily_0.02\Figures');
saveas(figure(1),[pwd '/Figures/mapSD' num2str(timeind) '.png']);        
saveas(figure(2),[pwd '/Figures/mapET' num2str(timeind) '.png']);        
saveas(figure(3),[pwd '/Figures/mapSrz' num2str(timeind) '.png']);        
saveas(figure(4),[pwd '/Figures/meascompSD_' num2str(row0) '_' num2str(col0) '.png']);        
saveas(figure(5),[pwd '/Figures/meascompSrz_' num2str(row0) '_' num2str(col0) '.png']);        
saveas(figure(6),[pwd '/Figures/meascompET_' num2str(row0) '_' num2str(col0) '.png']);
saveas(figure(7),[pwd '/Figures/compSD_' num2str(row) '_' num2str(col) '.png']);
saveas(figure(8),[pwd '/Figures/compSrz_' num2str(row) '_' num2str(col) '.png']);
saveas(figure(9),[pwd '/Figures/compET_' num2str(row) '_' num2str(col) '.png']);
saveas(figure(10),[pwd '/Figures/compQ.png']);
pause
close(figure(1));
close(figure(2));
close(figure(3)); 
close(figure(4));
close(figure(5));
close(figure(6));
close(figure(7));
close(figure(8));
close(figure(9)); 
close(figure(10));
