% MOD_WET_wrapper.m (script)
% Written by Steve Margulis 
%
% This is a "wrapper" script that allows MOD-WET to be run for different
% inputs. This can be useful for running the model to perform sensitivity 
% tests or in an ensemble context for open loop simulations and/or adapted 
% for use with a data assimilation method.
%
% Note: The MOD-WET model can simply be run by executing the
% "MOD_WET_model_driver" function (without any input arguments). Here we 
% are just modifying inputs through the optional "restart_info" input
% variables.
%
% For simplicity all outputs and perturbed inputs are saved to file for
% later use.
%
% This code assumes that the nominal parameters, initial conditions, and
% meteorological forcing inputs are those specified in or linked to in the
% input files in the MOD-WET path pointed to below. For the specified path
% (i.e. "toolbox_path" below), these would include:
%
% initial conditions specified in:
% "toolbox_path"/chapter11/initialize_model.m
%
% parameters specified in:
% "toolbox_path"/chapter11/set_static_physical_parameters.m
%
% nominal meteorological input file specified in:
% "toolbox_path"/set_control_parameters.m
% 
% Based on these nominal inputs, they are perturbed as shown below before
% being input to a realization of MOD-WET. So for an ensemble of inputs,
% there will be an ensemble of outputs.
clear all;close all

%% Specify path to root directory for the MOD-WET functions
toolbox_path='C:/Users/BANZH/Downloads/2019spring/CEE steve/MOD_WET_2017a/MOD_WET_2017a/MOD_WET_src';
% Add path
addpath(genpath(toolbox_path))

% Root directory where to store outputs
output_root_dir='C:/Users/BANZH/Downloads/2019spring/CEE Steve/Project_outputs/open_loop_ensemble1/';
mkdir(output_root_dir);
%% --------- SPECIFY INPUTS HERE ----------------------------------------%%
% Note: Additional inputs are specified in: model_input_uncertainty_parameters.m 
% and in the measurement model: measurement_model.m

% Specify whether to invoke some plots (set to 0 to skip plotting below
% or to 1 to perform plotting)
plot_flag=0;
% Number of ensemble members to run (set to 1 if running the true simulation)
N_ensemble_members=20; % generate open-loop ensemble

% Set random seeds for model input uncertainties (allows control so that
% the true is repeatable)
rng('default');
seed1=1000;
seed2=1001;
%% Define MOD-WET states
% These are all of the states that MOD-WET propagates and not necessarily
% all of those that will be updated via DA.
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
uncertain_state_numbers=[3]; % here Srz and SD from list above.
% Updating SWE only:
%uncertain_state_numbers=[5]; % here SWE from list above.

% Store names of states being treated as uncertain/updated
for ii=1:length(uncertain_state_numbers)
    uncertain_states(ii).names=mod_wet_states((uncertain_state_numbers(ii))).names;
end

% Define uncertain time-invariant parameters by name:
% These are not updated, but add uncertainty to the states via the model
% physics
uncertain_params(1).names='T0'; 
uncertain_params(2).names='K0'; % these two parameters are really coupled, but can control this by making them perfectly correlated below
%
% Define uncertain time-varying forcing by name:
% Meteorological inputs that are treated as uncertain. Note: As currently
% constructed MOD-WET takes in only a single set of forcings and then
% distributes the forcing spatially within the model. This makes it
% difficult to treat these inputs as random spatial fields. So instead they
% are just treated as correlated single perturbations that can vary in
% time, but are deterministically distributed in space according to the 
% MOD-WET disaggregation equations.
uncertain_forcings(1).names='PPT';
%% Call initialization codes to get nominal parameters for model simulation. 
% This first part is just a simple mechanism to load the nominal 
% parameters/initial from the existing MOD-WET parameter/initial condition 
% functions. The desired variable values should be set in those functions 
% and are loaded here. In other words, the nominal inputs should already be
% set in the normal MOD-WET inputs files. By running those functions you 
% will be storing all of the parameters into the "restart_info" structure
% array. The nominal values can then be perturbed and replaced as inputs 
% for a simulation (or ensemble of simulations) via the "restart_info" 
% variable that can be provided as an input to MOD-WET. Beyond
% parameter/forcing perturbations, this can also be used to set the
% start/end time of the simulation, etc. These will be useful in being able
% to use MOD-WET in data assimilation approaches.

% load nominal control parameters (i.e. start day, time step, number of 
% days, etc.)
[restart_info.control_params]=set_control_parameters;
% This is to make sure hourly maps are output (assumes forcing is provided
% at 0.25 hour interval)
restart_info.control_params.map_frq2store=4;
% load nominal static physical parameters (soil parameters, TOPMODEL 
% parameters, etc.
[restart_info.params]=set_static_physical_parameters(restart_info.control_params);
% pre-allocation of variables
[restart_info.states,~,~,~]= ...
                pre_allocate_variables(restart_info.control_params,restart_info.params);
% load nominal initial conditions for state variables (i.e. "old_vars")
[~,~,restart_info.old_vars]= ...
        initialize_model(restart_info.control_params,restart_info.params,restart_info.states);
% Do load nominal initial conditions in the loop, because the multiplier is
% not 1 mean, so it will have cumulative effects.
%% Grab domain-specific parameters
% Grab mask of pixel indices (the length of this vector is the size of 
% the state vector per variable.
pixel_indices=find(restart_info.params.static_maps.mask==1);
Npix=length(pixel_indices);
%% set some time series variable and map matrix
tsSD=nan(20,1465);
tsET=nan(20,1464);
tsQ=nan(20,1464);
tsSrz=nan(20,1465);
mapSD=nan(26,19,1464,20);
mapET=nan(26,19,1464,20);
mapSrz=nan(26,19,1464,20);
%% Set simulation control parameters:
% This might be useful for calling the model for different time
% periods (i.e. when using a filter). These are just copied from the
% set_control_parameters.m.
% Replace the starting day of the simulation (0 is the beg. of the WY)
restart_info.control_params.start_day=183; % start simulation
% Replace the length of the simulation and related input parameters.
restart_info.control_params.n_days=61; % number of days to simulate
restart_info.control_params.ntime=1./restart_info.control_params.dt*24.*restart_info.control_params.n_days;
% Length of time series to store (needs to be an integer)
restart_info.control_params.nt=restart_info.control_params.ntime/restart_info.control_params.timeseries_frq2store; 
% Set index to initialize forcings 
restart_info.control_params.start_time=1+...     
    restart_info.control_params.start_day*24./restart_info.control_params.dt;

%% Loop through an ensemble of simulations (with different inputs)        
for i_ensemble=1:N_ensemble_members

    %% Here the control parameters, physical parameters, initial conditions, and forcing can be perturbed:

    %% An example of changing the output file name (in order to store different files for different ensemble members
    % Replace output filename with one specific to this ensemble member
    restart_info.control_params.output_filename=[output_root_dir 'ensemble_member_' num2str(i_ensemble) '_output.mat'];
    %% Perturb. uncertain inputs:
    % load nominal forcings
    eval(['load ' '''' restart_info.control_params.met_data_filename ''''])

    %% Call function that generates perturbations based on specified statistics
    [IC_perturbations,params_perturbations,forcing_perturbations]=...    
    model_input_uncertainty_parameters(uncertain_states,uncertain_params,...
    uncertain_forcings,restart_info,length(PPT),seed1,seed2);
    % ---- Perturb initial conditions (based on statistics above):    
    % Loop through states that will be perturbed
    for ii=1:length(uncertain_states)
        if (plot_flag==1)
            figure;
            subplot(1,3,1);eval(['imagesc(restart_info.old_vars.' uncertain_states(ii).names '0)']);axis image;
            colorbar;title(['Nominal ' uncertain_states(ii).names])
            subplot(1,3,2);imagesc(IC_perturbations(:,:,ii));axis image; 
            colorbar;title(['Perturbation for ' uncertain_states(ii).names])
            subplot(1,3,3);eval(['imagesc(IC_perturbations(:,:,ii).* restart_info.old_vars.' uncertain_states(ii).names '0)'])
            axis image; colorbar;title(['Perturbed ' uncertain_states(ii).names])
            pause
        end
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

    % ---- Perturb time-invariant parameters (based on statistics above):
    % Loop through parameters that will be perturbed
    for ii=1:length(uncertain_params)
        if (plot_flag==1)
            figure;
            subplot(1,3,1);eval(['imagesc(restart_info.params.static_maps.' uncertain_params(ii).names ')']);axis image;
            colorbar;title(['Nominal ' uncertain_params(ii).names])
            subplot(1,3,2);imagesc(IC_perturbations(:,:,ii));axis image; 
            colorbar;title(['Perturbation for ' uncertain_params(ii).names])
            subplot(1,3,3);eval(['imagesc(params_perturbations(:,:,ii).*restart_info.params.static_maps.' uncertain_params(ii).names ')'])
            axis image; colorbar;title(['Perturbed ' uncertain_params(ii).names])
            pause
        end
        command=['restart_info.params.static_maps.' uncertain_params(ii).names ' = params_perturbations(:,:,ii).* restart_info.params.static_maps.' uncertain_params(ii).names ';'];
        eval(command)
    end
    % Specify the desired means and SIGMA of multiplicative perturbations
    %desired_mean = 1; % desired mean
    %desired_stddev= 0.2; % desired standard deviations (or coeffs. of variation if mean=1) 
    % Generate random multipliers based on the given mean and stddev
    %LN_factor=mvnrnd(desired_mean,desired_stddev);
    %factor=exp(LN_factor);
    % Perturb d_rz 
    %restart_info.params.d_rz=restart_info.params.d_rz*factor;
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
        if (plot_flag==1)
            figure;
            subplot(3,1,1)
            eval(['plot(' uncertain_forcings(ii).names ')']);
            ylabel(uncertain_forcings(ii).names);axis tight;grid
            subplot(3,1,2)
            eval(['plot(forcing_perturbations(ii,:))']);
            ylabel('Perturb.');axis tight;grid
            subplot(3,1,3)
            eval(['plot(' uncertain_forcings(ii).names '.*forcing_perturbations(ii,:))']);
            ylabel(['Perturb. ' uncertain_forcings(ii).names]);axis tight;grid
            pause
        end
        command=[uncertain_forcings(ii).names ' = ' uncertain_forcings(ii).names '.* forcing_perturbations(ii,:);' ];
        eval(command)
    end

    % Write forcings to file that will be loaded by MOD-WET
    restart_info.control_params.met_data_filename=[output_root_dir 'true_met_forcing_input.mat'];
    eval(['save ' '''' restart_info.control_params.met_data_filename '''' ' time gage_elev PPT SW Ta qa U Psfc'])
    
    %% Call MOD-WET with perturbed inputs
    MOD_WET_model_driver(restart_info)
    %% Load outputs and plot some states/fluxes
    eval(['load ' '''' restart_info.control_params.output_filename ''''])
    tsSD(i_ensemble,:)=states.time_series.SD;
    tsET(i_ensemble,:)=fluxes.time_series.ET;
    tsQ(i_ensemble,:)=fluxes.time_series.outlet_hydrograph;
    tsSrz(i_ensemble,:)=states.time_series.Srz;
    mapSD(:,:,:,i_ensemble)=states.maps.SD;
    mapET(:,:,:,i_ensemble)=fluxes.maps.ET;
    mapSrz(:,:,:,i_ensemble)=states.maps.Srz;
end
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/tsQ.mat','tsQ')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/tsET.mat','tsET')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/tsSD.mat','tsSD')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/mapET.mat','mapET')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/mapSD.mat','mapSD')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/mapSrz.mat','mapSrz')
save('C:/Users/BANZH/Downloads/2019spring/CEE steve/Project_outputs/open_loop_ensemble/tsSrz.mat','tsSrz')
DOY=183:1:243;
figure(1);hold on
pixSD=squeeze(mean(mapSD(14,17,:,:),4));
plot(DOY,mean(reshape(pixSD(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(1);hold on
    plot(DOY,mean(reshape(squeeze(mapSD(14,17,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('SD (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
figure(2);hold on
pixSD=squeeze(mean(mapSD(19,3,:,:),4));
plot(DOY,mean(reshape(pixSD(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(2);hold on
    plot(DOY,mean(reshape(squeeze(mapSD(19,3,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('SD (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')

figure(3);hold on
pixSD=squeeze(mean(mapSD(8,10,:,:),4));
plot(DOY,mean(reshape(pixSD(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(3);hold on
    plot(DOY,mean(reshape(squeeze(mapSD(8,10,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('SD (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
        
figure(4);hold on
pixSrz=squeeze(mean(mapSrz(14,17,:,:),4));
plot(DOY,mean(reshape(pixSrz(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(4);hold on
    plot(DOY,mean(reshape(squeeze(mapSrz(14,17,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('Srz (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
figure(5);hold on
pixSrz=squeeze(mean(mapSrz(19,3,:,:),4));
plot(DOY,mean(reshape(pixSrz(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(5);hold on
    plot(DOY,mean(reshape(squeeze(mapSrz(19,3,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('Srz (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')

figure(6);hold on
pixSrz=squeeze(mean(mapSrz(8,10,:,:),4));
plot(DOY,mean(reshape(pixSrz(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(6);hold on
    plot(DOY,mean(reshape(squeeze(mapSrz(8,10,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('Srz (m)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')

figure(7);hold on
pixET=squeeze(mean(mapET(14,17,:,:),4));
plot(DOY,mean(reshape(pixET(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(7);hold on
    plot(DOY,mean(reshape(squeeze(mapET(14,17,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('ET (m/h)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')
figure(8);hold on
pixET=squeeze(mean(mapET(19,3,:,:),4));
plot(DOY,mean(reshape(pixET(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(8);hold on
    plot(DOY,mean(reshape(squeeze(mapET(19,3,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('ET (m/h)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')

figure(9);hold on
pixET=squeeze(mean(mapET(8,10,:,:),4));
plot(DOY,mean(reshape(pixET(1:1464,1),24,61)),'k','LineWidth',3)
for i=1:20
    figure(9);hold on
    plot(DOY,mean(reshape(squeeze(mapET(8,10,:,i)),24,61)),'c');
end
        xlabel('DOWY');ylabel('ET (m/h)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','south')









































figure(4);hold on
meanET=mean(tsET,1);
plot(DOY,mean(reshape(meanET(1,1:1464),24,61)),'k','LineWidth',3)
for i=1:20
    figure(4);hold on
    plot(DOY,mean(reshape(tsET(i,1:1464),24,61)),'c')
end
        xlabel('DOWY');ylabel('ET (m/h)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','southeast')
        

figure(11);hold on
meanQ=mean(tsQ,1);
plot(DOY,mean(reshape(meanQ(1,1:1464),24,61)),'k','LineWidth',3)
for i=1:20
    figure(11);hold on
    plot(DOY,mean(reshape(tsQ(i,1:1464),24,61)),'c')
end
        xlabel('DOWY');ylabel('outlet hydrograph (m^3/s)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','southeast')
figure(11);hold on
meanSWE=mean(tsSWE,1);
plot(DOY,mean(reshape(meanSWE(1,1:1464),24,61)),'k','LineWidth',3)
for i=1:20
    figure(11);hold on
    plot(DOY,mean(reshape(tsSWE(i,1:1464),24,61)),'c')
end
        xlabel('DOWY');ylabel('outlet hydrograph (m^3/s)')
        set(gca,'FontSize',18,'YDir','normal');grid;
        legend('Ensemble mean','Each Ensemble','Location','southeast')
figure(4);hold on
stdSD=std(mapSD,0,4);
imagesc(stdSD(:,:,61));axis image;colorbar;
    caxis([min(min(stdSD(:,:,61))),max(max(stdSD(:,:,61)))])
    title(['ensemble std of SD at DOWY=243'])
    set(gca,'FontSize',18)  

figure(5);hold on
stdET=std(mapET,0,4);
imagesc(stdET(:,:,61));axis image;colorbar;
    caxis([min(min(stdET(:,:,61))),max(max(stdET(:,:,61)))])
    title(['ensemble std of LE at DOWY=243'])
    set(gca,'FontSize',18)  