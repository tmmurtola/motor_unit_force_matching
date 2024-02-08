classdef MUpool
    % MUpool is a motor unit pool object implementing the different MU pool
    % structures in
    %   Murtola & Richards (2024) ...

    properties
        tot_number
        dt
        name
        opts
    end

    methods
        function obj = MUpool(n,dt,opts)
            %MUPOOL(n,dt,opts) Construct an instance of this class
            %   Inputs:
            %       n (scalar)      number of motor units in pool
            %       dt (scalar)     size of timestep (in sec)
            %       opts (struct)   options for motor unit pool
            %   Output:
            %       MUpool object

            obj.tot_number = n;
            obj.dt = dt;
            obj.opts = opts;
        end

        function activation = drive2activation(obj,drive)
            % DRIVE2ACTIVATION(obj, drive) Compute normalised activation 
            % level of MUs in pool corresponding to neural drive signal
            %   Input:
            %       obj         MUpool object (with n MUs)
            %       drive       neural drive signal as a row vector (1 x k)
            %   Output:
            %       activation  matrix of MU activation signals (n x k)
            %                   (row i = activation for MU i, col j = time
            %                   step j)

            frs = drive2rate(obj,drive);
            us = rate2excitation(obj,frs);
            activation = excitation2activation(obj,us);

        end

        function firing_rates = drive2rate(obj,drive)
            % DRIVE2RATE(obj, drive) Compute firing rates of MUs
            % in pool corresponding to neural drive signal
            %   Input:
            %       obj         MUpool object (with n MUs)
            %       drive       neural drive signal as a row vector (1 x k)
            %   Output:
            %       firing_rates  matrix of firing rates (n x k)
            %                   (row i = firing rate for MU i, col j = time
            %                   step j)

            switch obj.opts.rate_function

                case 'linear' % linear, capped rate function
                    m = (obj.opts.rate_function_params.max_rate - obj.opts.rate_function_params.min_rate)./(obj.opts.rate_function_params.max_threshold-obj.opts.recr_threshold);
                    firing_rates = obj.opts.rate_function_params.min_rate + m.*(drive-obj.opts.recr_threshold);
                    % remove rates below recruitment threshold
                    firing_rates(drive < obj.opts.recr_threshold) = 0;
                    % cap at maximum firing rate
                    if length(obj.opts.rate_function_params.max_rate) == 1  % apply single max_rate to all MUs
                        firing_rates(drive >= obj.opts.rate_function_params.max_threshold) = obj.opts.rate_function_params.max_rate;
                    else   % each MU has its own max_rate
                        temp_max = repmat(obj.opts.rate_function_params.max_rate,1,size(firing_rates,2));
                        firing_rates(drive >= obj.opts.rate_function_params.max_threshold) = temp_max(drive >= obj.opts.rate_function_params.max_threshold);
                    end

                case 'log' % logarithm based rate function
                    firing_rates = 1./(1./obj.opts.rate_function_params.max_rate-obj.opts.rate_function_params.ref_drive./obj.opts.recr_threshold.*log(1-obj.opts.recr_threshold./drive));
                    % remove rates below recruitment threshold
                    firing_rates(drive<=obj.opts.recr_threshold) = 0;

                otherwise
                    error('Unknown rate function. Check obj.opts.rate_function.')

            end
        end

        function excitations = rate2excitation(obj,frs)
            % RATE2EXCITATION(obj, frs) Compute excitation impulse trains 
            % for MUs in pool corresponding to firing rate signals
            %   Input:
            %       obj         MUpool object (with n MUs)
            %       frs         firing rate matrix (n x k)
            %                   (row i = firing rate for MU i, col j = time
            %                   step j)
            %   Output:
            %       excitations  matrix of excitation signals (n x k)
            %                   (row i = activation for MU i, col j = time
            %                   step j)

            excitations = zeros(size(frs));
            
            % compute firing times one MU at a time
            for i = 1:obj.tot_number

                % find time window when firing rate > 0 
                start_ind = find(frs(i,:)>0,1,'first');
                end_ind = find(frs(i,:)>0,1,'last');
                if end_ind == start_ind % frs never goes back to zero
                    end_ind = length(frs(i,:));
                end

                % desired inter-impulse time
                IITs = 1./frs(i,:);

                % first spike when frs>0 for the first time
                excitations(i,start_ind) = 1;
                time_since_fired = obj.dt;
                
                % loop through time steps of interest
                for j = (start_ind+1):end_ind
                    if IITs(j)<= time_since_fired  % fire at current time step
                        excitations(i,j) = 1;  
                        time_since_fired = time_since_fired - IITs(j) + obj.dt; % add discrete time correction
                    else    % no fire, just update time_since_fired
                        time_since_fired = time_since_fired + obj.dt;
                    end
                end
            end

        end

        function activation = excitation2activation(obj,u)
            % EXCITATION2ACTIVATION(obj, u) Compute normalised activation 
            % level of MUs in pool corresponding to excitation 
            % impulse trains
            %   Input:
            %       obj         MUpool object (with n MUs)
            %       u           excitation signal matrix (n x k)
            %                   (row i = excitation for MU i, col j = time
            %                   step j)
            %   Output:
            %       activation  matrix of MU activation signals (n x k)
            %                   (row i = activation for MU i, col j = time
            %                   step j)

            activation = zeros(size(u));
            % intermediate states for activation dynamics
            a1 = zeros(obj.tot_number,1);
            a2 = zeros(obj.tot_number,1);
            a3 = zeros(obj.tot_number,1);

            for i = 1:size(u,2)
                a1 = a1 + obj.dt./obj.opts.alpha(:,1).*(u(:,i)-(obj.opts.beta(:,1)+(1-obj.opts.beta(:,1)).*u(:,i)).*a1);
                a2 = a2 + obj.dt./obj.opts.alpha(:,2).*(a1    -(obj.opts.beta(:,2)+(1-obj.opts.beta(:,2)).*a1)    .*a2);
                a3 = a3 + obj.dt./obj.opts.alpha(:,3).*(a2    -(obj.opts.beta(:,3)+(1-obj.opts.beta(:,3)).*a2)    .*a3);
                activation(:,i) = a3;
            end
            activation = activation./obj.opts.a_max;
        end
    end
end