classdef MUpool
    % MUpool is a motor unit pool object implementing the different MU pool
    % structures in
    %   Murtola & Richards (2024) Matching dynamically varying forces with multi-motor-unit muscle models: A simulation study
    %
    % (c) 2024 Tiina Murtola/RVC

    properties
        tot_number
        dt
        name
        opts
        a1
        a2
        a3
        time_since_fired
        last_IIT
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
            obj.a1 = zeros(n,1);
            obj.a2 = zeros(n,1);
            obj.a3 = zeros(n,1);
            obj.time_since_fired = inf(n,1);
            obj.last_IIT = inf(n,1);
        end

        function obj = reset(obj)
            %RESET(obj) Reset state variables in MU pool
            %   Input:
            %       obj         MUpool object (with n MUs)
            %   Output:
            %       obj         MUpool object (with n MUs)
            obj.a1 = zeros(obj.tot_number,1);
            obj.a2 = zeros(obj.tot_number,1);
            obj.a3 = zeros(obj.tot_number,1);
            obj.time_since_fired = inf(obj.tot_number,1);
            obj.last_IIT = inf(obj.tot_number,1);
        end

        function activation = drive2activation(obj,drive)
            % DRIVE2ACTIVATION(obj, drive) Compute normalised activation 
            % level of MUs in pool corresponding to neural drive signal.
            % Ignores intial state of MUs.
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

        function [activation,obj] = drive2activationStep(obj,drive)
            % DRIVE2ACTIVATIONSTEP(obj, drive) Compute normalised activation 
            % level of MUs in pool corresponding to neural drive signal,
            % starting from current state of MUs and returning updated MU
            % state
            %   Input:
            %       obj         MUpool object (with n MUs)
            %       drive       neural drive signal as a row vector (1 x k)
            %   Output:
            %       activation  matrix of MU activation signals (n x k)
            %                   (row i = activation for MU i, col j = time
            %                   step j)
            %       obj         MUpool object (with n MUs)

            frs = drive2rate(obj,drive);
            [us,obj] = rate2excitationStep(obj,frs);
            [activation,obj] = excitation2activationStep(obj,us);

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

                % find time window when firing rate >= min_rate 
                start_ind = find(frs(i,:)>0,1,'first');
                end_ind = find(frs(i,:)>0,1,'last');
                if end_ind == start_ind % frs never goes back to zero
                    end_ind = length(frs(i,:));
                end

                % desired inter-impulse times
                IITs = 1./frs(i,:);

                % first spike when frs>=min_rate for the first time
                excitations(i,start_ind) = 1;
                time_since_fired = obj.dt;
                
                % loop through time steps of interest
                for j = (start_ind+1):end_ind
                    if IITs(j)<= time_since_fired  % fire at current time step
                        excitations(i,j) = 1;
                        if isinf(IITs(j-1))    % first & derecruited
                            time_since_fired = 0;
                        else
                            time_since_fired = obj.dt - (IITs(j-1)-time_since_fired+obj.dt)/(1-(IITs(j)-IITs(j-1))/obj.dt); % add discrete time correction
                        end
                    end    
                    % update time_since_fired
                    time_since_fired = time_since_fired + obj.dt;                    
                end
            end

        end

        function [excitations, obj] = rate2excitationStep(obj,frs)
            % RATE2EXCITATIONSTEP(obj, frs) Compute excitation impulse trains 
            % for MUs in pool corresponding to firing rate signals,             
            % starting from current state of MUs and returning updated MU
            % state
            %   Input:
            %       obj         MUpool object (with n MUs)
            %       frs         firing rate matrix (n x k)
            %                   (row i = firing rate for MU i, col j = time
            %                   step j)
            %   Output:
            %       excitations  matrix of excitation signals (n x k)
            %                   (row i = activation for MU i, col j = time
            %                   step j)
            %       obj         MUpool object (with n MUs)

            excitations = zeros(size(frs));
            
            % compute firing times one MU at a time
            for i = 1:obj.tot_number

                % desired inter-impulse times
                IITs = 1./frs(i,:);
                
                % loop through time steps
                for j = 1:size(frs,2)
                    if ~isinf(IITs(j))
                        if IITs(j)<= obj.time_since_fired(i)  % fire at current time step
                            excitations(i,j) = 1;
                            if isinf(obj.last_IIT(i))    % first & derecruited
                                obj.time_since_fired(i) = 0;
                            else
                                obj.time_since_fired(i) = obj.dt - (obj.last_IIT(i)-obj.time_since_fired(i)+obj.dt)/(1-(IITs(j)-obj.last_IIT(i))/obj.dt); % add discrete time correction
                            end
                        end
                    end
                    % update time_since_fired & last IIT
                    obj.time_since_fired(i) = obj.time_since_fired(i) + obj.dt;
                    obj.last_IIT(i) = IITs(j);
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

        function [activation,obj] = excitation2activationStep(obj,u)
            % EXCITATION2ACTIVATION(obj, u) Compute normalised activation 
            % level of MUs in pool corresponding to excitation 
            % impulse trains starting from current state of MUs and 
            % returning updated MU state
            %   Input:
            %       obj         MUpool object (with n MUs)
            %       u           excitation signal matrix (n x k)
            %                   (row i = excitation for MU i, col j = time
            %                   step j)
            %   Output:
            %       activation  matrix of MU activation signals (n x k)
            %                   (row i = activation for MU i, col j = time
            %                   step j)
            %       obj         MUpool object (with n MUs)

            activation = zeros(size(u));

            for i = 1:size(u,2)
                obj.a1 = obj.a1 + obj.dt./obj.opts.alpha(:,1).*(u(:,i)  -(obj.opts.beta(:,1)+(1-obj.opts.beta(:,1)).*u(:,i)).*obj.a1);
                obj.a2 = obj.a2 + obj.dt./obj.opts.alpha(:,2).*(obj.a1  -(obj.opts.beta(:,2)+(1-obj.opts.beta(:,2)).*obj.a1).*obj.a2);
                obj.a3 = obj.a3 + obj.dt./obj.opts.alpha(:,3).*(obj.a2  -(obj.opts.beta(:,3)+(1-obj.opts.beta(:,3)).*obj.a2).*obj.a3);
                activation(:,i) = obj.a3;
            end
            activation = activation./obj.opts.a_max;
        end
    end
end