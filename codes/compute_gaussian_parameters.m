function [gaussian_mu,gaussian_covariance_matrix]=...
    compute_gaussian_parameters(desired_mean_vector,...
    desired_stddev_vector,desired_corr_matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function transform lognormal inputs into Gaussian parameters by 
% passing the desired mean, standard deviations, and correlation matrix 
% for the lognormal random vector and returning the transformed (Gaussian) 
% mu parameter and covariance matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of variables in vector
nvar=length(desired_mean_vector);

% Compute log-normal parameters
desired_variances = desired_stddev_vector.^2; % variance  
gaussian_mu = log((desired_mean_vector.^2) ./ sqrt(desired_variances + ...
    desired_mean_vector.^2)); % parameter of the associated normal dist.
% parameter of the associated normal dist.
sigma = sqrt(log(desired_variances ./ (desired_mean_vector.^2)+1));

% Transform to covariance matrix for normal dist.
gaussian_covariance_matrix=NaN(nvar,nvar);
for i=1:nvar
    for j=1:nvar
        gaussian_covariance_matrix(i,j)=log( desired_corr_matrix(i,j)* ... 
            sqrt( exp(sigma(i)^2)-1)*sqrt(exp(sigma(j)^2)-1) + 1 );        
    end
end

return