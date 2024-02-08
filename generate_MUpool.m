function mup = generate_MUpool(n,dt,Fmax,Rf,alpha,beta,threshold_mode,threshold_params,rate_func,rate_params)
% GENERATE_MUPOOL Initialises a MUpool object from given parameters

opts = struct();

mui = (1:n)'; % MU index

% -- MU strenghts -- %

if n == 1
    fmax = Fmax;
else
    fmax = Rf.^((0:n-1)'/(n-1));
    fmax = fmax./sum(fmax) * Fmax;
end

opts.max_force = fmax;


% -- Rate function -- %

switch rate_func

    case 'linear_const_max' % linear + constant r_max, rate_params = [min_rate max_rate]
        opts.rate_function = 'linear';
        opts.rate_function_params.min_rate = min(rate_params);
        opts.rate_function_params.max_rate = max(rate_params);
    case 'linear_lindecr' % linear +  [min_rate max_rate]
        opts.rate_function = 'linear';
        opts.rate_function_params.min_rate = min(rate_params);
        opts.rate_function_params.max_rate = max(rate_params).*(1-(mui-1)./(n-1)*0.25);
    case 'linear_expdecr' % linear + exp decrease in r_max, rate_params [min_rate max_rate]
        opts.rate_function = 'linear';
        opts.rate_function_params.min_rate = min(rate_params);
        opts.rate_function_params.max_rate = max(rate_params).*mui.^-0.05;

    case 'log_const_max' % logarithmic + constant r_max, rate_params = max_rate
        opts.rate_function = 'log';
        opts.rate_function_params.ref_drive = 0.1*2.4.^((mui./n).^1.47);
        opts.rate_function_params.max_rate = rate_params(1);
    case 'log_lindecr' % logarithmic + lin decrease in r_max, rate_params = max_rate
        opts.rate_function = 'log';
        opts.rate_function_params.ref_drive = 0.1*2.4.^((mui./n).^1.47);
        opts.rate_function_params.max_rate = rate_params(1).*(1-(mui-1)./(n-1).*0.25);
    case 'log_expdecr' % logarithmic + exp decrease in r_max, rate_params = max_rate
        opts.rate_function = 'log';
        opts.rate_function_params.ref_drive = 0.1*2.4.^((mui./n).^1.47);
        opts.rate_function_params.max_rate = rate_params(1).*mui.^-0.05;

    otherwise
        error('Unknown rate function.')

end


% -- Activation dynamics parameters -- %

if size(alpha,1) == n
    opts.alpha = alpha;
elseif size(alpha,1) >1
    error('Number of rows of alpha must be 1 or n.')
else
    opts.alpha = repmat(alpha,n,1);
end

if size(beta,1) == n
    opts.beta = beta;
elseif size(beta,1) >1
    error('Number of rows of alpha must be 1 or n.')
else
    opts.beta = repmat(beta,n,1);
end

beta_tot = prod(opts.beta,2);
a_max = 1./(beta_tot./(opts.rate_function_params.max_rate.*dt) + 1 -beta_tot);
opts.a_max = a_max;

% -- Recruitment thresholds -- %

switch threshold_mode
    case 'exponential'
        first_recruitment = threshold_params(1);
        last_recruitment = threshold_params(2);
        recruitment_bw = Fmax-last_recruitment;
        RR = last_recruitment/first_recruitment;

        min_thr = first_recruitment * RR.^((0:n-1)'/(n-1));
        max_thr = min_thr+recruitment_bw;

    case 'matched' 
        first_recruitment = threshold_params(1);
        overlap = 1 - (threshold_params(2)-threshold_params(1))/(Fmax-fmax(end));
        recruitment_bw = fmax;

        min_thr = zeros(n,1);
        min_thr(1) = first_recruitment;
        for i = 2:n
            min_thr(i)= min_thr(i-1)+recruitment_bw(i-1)*(1-overlap);
        end
        max_thr = min_thr+recruitment_bw;

    case 'mixed_exponential'
        c1 = threshold_params(1);
        c2 = 29.1951;
        c3 = 1.8333;
        muir = (0:n-1)'./(n-1);
        Df = (threshold_params(2)-c2)/c1;
        min_thr = c1* Df.^(muir.^c3) + c2*muir;

end

opts.recr_threshold = min_thr;

if exist("max_thr","var")
    opts.rate_function_params.max_threshold = max_thr;
end

% -- generate MUpool object -- %

mup = MUpool(n,dt,opts);

end