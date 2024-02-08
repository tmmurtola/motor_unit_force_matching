function [npools, mups, Fmax] = muscle_MU_settings(muscle, dt)

% MU pool parameters by muscle

switch muscle

    case 'TA'
        Fmax = 400;     % desired Fmax
        n = 400;        % number of MUs

        min_rate = 0;   % minimum firing rate (for all MUs)
        max_rate = 25;  % maximum firing rate (MU 1)
        Rf = 100;       % range of MU strengths

        alpha = [0.014 0.0356 0.0223];
        beta = [0.6103 0.7815 0.6871];

        npools = 9;
        thresholds = {'exponential','exponential','exponential','matched','matched','matched','mixed_exponential','mixed_exponential','mixed_exponential'};
        threshold_params = repmat([0.5 0.5*Fmax],npools,1);
        rate_funcs = {'linear_const_max','linear_lindecr','linear_expdecr','linear_const_max','linear_lindecr','linear_expdecr','log_const_max','log_lindecr','log_expdecr'};
        rate_params = [repmat([min_rate max_rate],6,1); repmat([max_rate nan],3,1)];
        pool_names ={'exp-lin const','exp-lin lin','exp-lin exp','matched-lin const','matched-lin lin','matched-lin exp','mixed-log const','mixed-log lin','mixed-log exp'};

    case 'SF'
        Fmax = 200;
        n = 200;

        min_rate = 0;
        max_rate = 50;
        Rf = 10;

        alpha = [0.018 0.04 0.025];
        beta = [0.6 0.8 0.7];

        npools = 9;
        thresholds = {'exponential','exponential','exponential','matched','matched','matched','mixed_exponential','mixed_exponential','mixed_exponential'};
        threshold_params = repmat([0.5 0.95*Fmax],npools,1);
        rate_funcs = {'linear_const_max','linear_lindecr','linear_expdecr','linear_const_max','linear_lindecr','linear_expdecr','log_const_max','log_lindecr','log_expdecr'};
        rate_params = [repmat([min_rate max_rate],6,1); repmat([max_rate nan],3,1)];
        pool_names ={'exp-lin const','exp-lin lin','exp-lin exp','matched-lin const','matched-lin lin','matched-lin exp','mixed-log const','mixed-log lin','mixed-log exp'};
end

% generate MU pools

for k = 1:npools
    mups(k)= generate_MUpool(n,dt,Fmax,Rf,alpha,beta,thresholds{k},threshold_params(k,:),rate_funcs{k},rate_params(k,:));
    mups(k).name = pool_names{k};
end
end