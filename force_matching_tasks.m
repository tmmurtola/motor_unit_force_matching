function [taskTargets,dt,normFlag] = force_matching_tasks(taskType,varargin)
% FORCE_MATCHING_TASKS generates force and FVL profiles for MU force
% matching tasks
%   Inputs:
%       taskType    'simulated_reach', 'trapezoidal', or 'elementary'
%       taskDur     (opt) task duration for elementary tasks (scalar)
%       lengthChange (opt) vector of relative muscle length changes for 
%                    elementary tasks
%       lengthOffset (opt) initial normalised muscle length (scalar)


tau_init = 50;      % length of initial pause before profiles start

switch taskType
    case 'simulated_reach'
        % load data structure
        load("reaching_data.mat","S");
        
        % extract time step (same for all simulations)
        dt = S(1).time(2)-S(1).time(1);

        cut_len = 1050;             % number of samples to include
        muscle_id = 1;              % muscle order in data: shoulder flx, shoulder ext, elbow flx, elbow ext, wrist flx, wrist ext

        taskTargets = nan(size(S,2),tau_init+cut_len,2);
        
        % extract force and FVL gains data
        for i = 1:size(S,2)
            taskTargets(i,:,1) = [zeros(1,tau_init) -S(i).frc(1:cut_len,muscle_id)'];
            taskTargets(i,:,2) = [ones(1,tau_init) S(i).glen(1:cut_len,muscle_id)'.*S(i).gvel(1:cut_len,muscle_id)'];
        end

        % fix issue with glen
        taskTargets(:,tau_init+1,2) = 1;

        normFlag = 0;               % force does not need scaling

    case 'trapezoidal'
        dt = 0.002;
        time = 0:dt:37;
        init_time = tau_init*dt;

        taskTargets = zeros(2,size(time,2));
        ramps = [6 10];             % duration of ramps for each task
        holds = [20 15];            % duration of holds for each task

        % compute force profiles
        for i = 1:2
            ramp_dur = ramps(i);
            hold_dur = holds(i);
            taskTargets(i,time>=init_time & time<(init_time+ramp_dur)) = 1/ramp_dur* (time(time>=init_time & time<(init_time+ramp_dur))-init_time);
            taskTargets(i,time>=(init_time+ramp_dur) & time<(init_time+ramp_dur+hold_dur)) = 1;
            taskTargets(i,time>=(init_time+ramp_dur+hold_dur) & time<(init_time+2*ramp_dur+hold_dur)) = 1-1/ramp_dur* (time(time>=(init_time+ramp_dur+hold_dur) & time<(init_time+2*ramp_dur+hold_dur))-(init_time+ramp_dur+hold_dur));
        end

        taskTargets(:,:,2) = 1;     % isometric: FVL = 1

        normFlag = 1;               % force is normalised

    case 'elementary'
        dt = 0.002;
        taskDur = varargin{1};
        lengthChange = varargin{2};
        lengthOffset = varargin{3};
        time = 0:dt:(taskDur+0.5);      % add extra time  

        m = length(lengthChange);  % number of tasks

        taskTargets = zeros(m,length(time),2);
        taskTargets(:,:,2) = 1;                         % temporary initial value
 
        for j = 1:length(lengthChange)
            Lplan = lengthOffset + [0 lengthChange(j)]; % beginning and end lengths
            [norm_frc, fvl_gain] = make_nonisom_elementary(taskDur,Lplan,dt,tau_init);
            taskTargets(j,1:length(norm_frc),1) = norm_frc;
            taskTargets(j,1:length(norm_frc),2) = fvl_gain;
            taskTargets(j,length(norm_frc)+1:end,2) = fvl_gain(end);
        end
        

        normFlag = 1;

    otherwise
        error('Unknown task type (input variable 1)')
end

end


function [force, fvl] = make_nonisom_elementary(Tplan, Lplan, dt, tau_init)

time = 0:dt:Tplan;
tau = time./Tplan;

w = 6*tau.^5-15*tau.^4+10*tau.^3;
dw = 1/Tplan*(6*5*tau.^4-15*4*tau.^3+10*3*tau.^2);
ddw = 1/Tplan/Tplan*(6*5*4*tau.^3-15*4*3*tau.^2+10*3*2*tau.^1);


muscle_len = [Lplan(1)*ones(1,tau_init) Lplan(1)+(Lplan(2)-Lplan(1))*w Lplan(2)*ones(1,tau_init)];
muscle_vel = [zeros(1,tau_init) (Lplan(1)-Lplan(2))*dw zeros(1,tau_init)];

acc_av = (ddw(1:floor(length(ddw)/2))-ddw(ceil(length(ddw)/2)+1:end))/2;

if Lplan(1)-Lplan(2)>=0 % contraction or isometric: muscle = agonist
    muscle_acc = [zeros(1,tau_init) acc_av zeros(1,length(muscle_vel)-length(acc_av)-tau_init)]./max(acc_av);
else
    muscle_acc = [zeros(1,length(muscle_vel)-length(acc_av)-tau_init) acc_av zeros(1,tau_init)]./max(acc_av);
end

force = muscle_acc;


b = [1.5 -1.1 0.53];    % active FL parameters
fl = @(x) exp(-abs((x.^b(2)-1)/b(3)).^b(1));

c = [4 1.8 30.24];  % FV parameters
fvc = @(x) (1-x)./(1+c(1)*x);   % contraction (positive x)
fve = @(x) c(2) - (c(2)-1)*(1+x)./(1-c(3)*x);   % elongantion (negative x)
fv = @(x) (x>0).*fvc(abs(x)) + (x<=0).*fve(-abs(x));

vmax = 1.6;

fvl = fl(muscle_len).*fv(muscle_vel/vmax);

end