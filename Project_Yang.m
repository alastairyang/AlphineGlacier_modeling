%% TITLE: The impact of changing climate at high elevation on Guoqu Glacier
% Glaciology Mid-term research project
%% Author: Donglai Yang
% Date: April 12, 2022

%% To run or not to run...?
% If RunModel = 1, you are going to sit and wait for probably 20 mins
% (equivalent to 2 youtube videos) since initializing the glacier to be of the right thickness takes 1. some
% trials and 2. a long time given how little the precipitation there is at
% ~5500 meter.
% If RunModel = 0, you are just loading the data I saved.
RunModel = 1;

%% Data
guoqu_table = readtable('GuoquGlacier_FlowlineData.csv');

dist = guoqu_table.distance;
H = flip(guoqu_table.THICKNESS_);
S = flip(guoqu_table.SurfaceDEM);
H(isnan(H)) = 0;
b = S - H;
b = movmean(b, 4);

termi_i = find(H==0, 1);
S(termi_i:end) = b(termi_i:end);

% Interpolate to a regular grid
dx = 10; % meter
x = dist(1):dx:dist(end);
zb = interp1(dist, b, x);
H  = interp1(dist, H, x);
% aggregate into a structure
geometry.dx = dx;
geometry.x = x;
geometry.zb = zb;

% Current length
L = dist(1) - dist(termi_i); % meter

% plot
figure(1);
% Define plot colors
rockcolor = [0.7 0.5 0.2];   % in RGB space; it's light brown
icecolor =  [0.7 0.9 1];     % in RGB space; it's light blue
ax(1) = subplot(2,1,1,'replace'); hold on
area(x/1e3,(zb+H)','facecolor',icecolor,'basevalue',5000);
area(x/1e3,zb,'facecolor',rockcolor,'basevalue',5000)
title(sprintf('Guoqu glacier with ELA at %1.0f m, year %1.0f',5570,1982))
xlabel('Horizontal Distance (km)','FontName','Times','FontSize',13)
ylabel('Elevation (m)','FontName','Times','FontSize',13)
print(gcf,'results/glacier.png','-dpng','-r300')

%% Parametric sweep
dbdzs = 0.00005:0.00004:0.0003; % m/yr/m
ELAs  = 5570;
N_dbdzs = length(dbdzs);
N_ELAs  = length(ELAs);
MSE = zeros(N_dbdzs, N_ELAs);
model_H = zeros(N_dbdzs*N_ELAs, length(zb));
parameters = zeros(N_dbdzs*N_ELAs, 2);

if RunModel == 1
    % for loop
    for i = 1:N_dbdzs
        for j = 1:N_ELAs
            parameter.dbdz = dbdzs(i);
            parameter.ELA  = ELAs(j);
            model_out = FlowLineM_init(geometry, parameter);
            MSE(i,j) = mean((model_out.H - H).^2);
            
            % get linear index
            i_ind = sub2ind([N_dbdzs,N_ELAs],i,j);
            model_H(i_ind,:) = model_out.H;
            
            % Parameters used
            parameters(i_ind,:) = [parameter.dbdz, parameter.ELA];
        end
    end
    save results/Initialize_Hs.mat model_H
    save results/parameters.mat parameters
    save results/MSE.mat MSE
else
    load results/Initialize_Hs.mat
    load results/parameters.mat
end
%% plots
% Plot thickness profiles
% Get output from the initalization
[~, min_i] = min(MSE);

figure('Position',[100,100,800,600]);
subplot(1,2,1)
plot(x, H,'-.','LineWidth',2); hold on
plot(x, model_H)
dbdz_strs = "db/dz = "+string(dbdzs);
legend([{'Observations'}, convertStringsToChars(dbdz_strs)])
title('Observed and simulated ice thickness','FontName','Times','FontSize',13)
xlabel('Along flow distance (m)','FontName','Times','FontSize',13)
ylabel('Thickness (m)','FontName','Times','FontSize',13)

% Plot MSE
subplot(1,2,2)
plot(parameters(:,1), MSE, '-*','LineWidth', 3)
xticks(parameters(:,1))
xlabel('db/dz','FontName','Times','FontSize',13)
ylabel('MSE','FontName','Times','FontSize',13)
title('Mean Squared error of ice thickness','FontName','Times','FontSize',13)

print(gcf, 'results/thickness_MSE.png','-dpng','-r300')
%% Climate perturbation: ELA moves up from 1986 - 2006
ddbdzdt = 0.000004:-0.0000002:-0.000008; % mass balance gradient does not change
dELAdt = 12:0.5:40; % m/yr, positive -> going up in altitude
% optimal dbdz_0 and ELA we found from initialization
ELA_init = 5570;
dbdz_init = parameters(3,1);
parameter.ELA = ELA_init;
parameter.dbdz = dbdz_init;
geometry.H_init = model_H(3,:);

% find out the range of accumulation
% we know 1982 record exists, 1986 does not
% so the deficit can be any year in between
dt_accum = 1986 - 1982;
b_accum = (dbdz_init*(zb + geometry.H_init - ELA_init))*dt_accum;
most_deficit = interp1(zb+geometry.H_init, b_accum, 5750);
least_deficit = 0; % happen to just melt 1986 layer

% model time: 0 to 20 years
% corresponding to 1986 - 2006
time = 20;

% run the model
H_perturb = zeros(length(dELAdt), length(ddbdzdt), length(x));
L_perturb = zeros(length(dELAdt), length(ddbdzdt));
accum_perturb = zeros(length(dELAdt), length(ddbdzdt));
melt_perturb = zeros(length(dELAdt), length(ddbdzdt));

if RunModel == 1
    count = 0;
    for i = 1:length(dELAdt)
        for j = 1:length(ddbdzdt)
            parameter.dELAdt  = dELAdt(i);
            parameter.ddbdzdt = ddbdzdt(j);
            pmodel_out = FlowLineM_perturb(geometry, parameter, time);
            H_perturb(i,j,:) = pmodel_out.H;
            L_perturb(i,j) = pmodel_out.L;
            accum_perturb(i,j) = pmodel_out.accum;
            melt_perturb(i,j) = pmodel_out.melt*dx;
            count = count + 1;
            disp(count)
        end
    end
    
    % Make ddbdzdt and dELAdt a meshgrid
    % for saving results
    [ddbdzdt_X, dELAdt_Y] = meshgrid(ddbdzdt, dELAdt);
    results.dELAdt = dELAdt_X;
    results.ddbdzdt = ddbdzdt_Y;
    results.H_pt = H_perturb;
    results.L_pt = L_perturb;
    results.accum_pt = accum_perturb;
    results.melt_pt  = melt_perturb;

    save results/perturb_results.mat results
else
    load results/perturb_results.mat
    % unpack the data
    dELAdt_X  = results.dELAdt;
    ddbdzdt_Y = results.ddbdzdt;
    H_perturb    = results.H_pt;
    L_perturb    = results.L_pt;
    accum_perturb = results.accum_pt;
    melt_perturb  = results.melt_pt;
end

%% Make plots

% get accum range
accum_range = accum_perturb < 0 & accum_perturb > -most_deficit;
new_accum = NaN(size(accum_perturb));
new_accum(accum_range) = accum_perturb(accum_range);
figure('Position',[100,100,700,700])
heatmap(ddbdzdt, dELAdt, new_accum)
xlabel('db/dz change rate ((m/s/m)/s)')
ylabel('ELA migration rate (m/s)')
print(gcf,'results/net_smb.png','-dpng','-r300')

% melt volume
melt_map = NaN(size(melt_perturb));
melt_map(accum_range) = melt_perturb(accum_range);

figure('Position',[100,100,700,700])
heatmap(ddbdzdt, dELAdt, melt_map);
xlabel('db/dz change rate ((m/s/m)/s)')
ylabel('ELA migration rate (m/s)')
print(gcf,'results/melt_19802006.png','-dpng','-r300')


%% Future projection
time = 120; % 120 year to 2100
N_MC = 100;

time_standard = 1:time;
% scenario 2: no change
if RunModel == 1 
    parameter.ddbdzdt = ddbdzdt(21);
    std_dELAdt2 = abs(dELAdt(22) - dELAdt(15))/4; % assuming the limits span 4 std
    mean_dELAdt2 = dELAdt(19);
    MC2 = randn(1,N_MC)*std_dELAdt2 + mean_dELAdt2;
    for i = 1:N_MC
        parameter.dELAdt = MC2(i);
        pj_model_out = FlowLineM_proj(geometry, parameter, time);
        results.time = pj_model_out.time_ax;
        results.melt = cumsum(pj_model_out.dmelt);
        results.Hleft = results.melt/sum(geometry.H_init*dx);
        results.L    = pj_model_out.L;
        MC2_results.(['No',num2str(i)]) = results;
    end
    
    % Interpolate to new time axis
    trialnames = fieldnames(MC2_results);
    MC_save = zeros(N_MC, length(time_standard));
    Hleft_save = zeros(N_MC, length(time_standard));
    Ls = zeros(N_MC);
    for i = 1:N_MC
        this_trial = MC2_results.(trialnames{i});
        melt_interp = interp1(this_trial.time, this_trial.melt, time_standard);
        Hleft = interp1(this_trial.time, this_trial.Hleft, time_standard);
        Ls(i) = this_trial.L;
        
        % save to a matrix for smaller file size
        MC_save(i,:) = melt_interp;
        Hleft_save(i,:) = Hleft;
    end
    % save
    save results/MC_save.mat MC_save
    save results/Hleft.mat Hleft
    save results/Ls.mat Ls

else
    load results/MC_save.mat
    load results/Hleft.mat
    load results/Ls.mat
end
%% Make a plot
% Assuming the glacier width is 1.5 km, which is the width of Guoqu Glacier
% in the middle of its trunk
width = 1500; % m
MeltIce_V = MC_save*width; % m^2


figure;
subplot(1,3,1)
plot(time_standard, MeltIce_V/1000^3)
xlabel('Time (year) from 1980','FontName','Times','FontSize',13)
ylabel('Melt ice volume (km^3)','FontName','Times','FontSize',13)
subplot(1,3,2)
plot(time_standard, Hleft_save)
xlabel('Time (year) from 1980','FontName','Times','FontSize',13)
ylabel('Percent of original ice volume','FontName','Times','FontSize',13)
subplot(1,3,3)
histogram(Ls(Ls>0))
xlabel('Glacier length from ice divide (m)','FontName','Times','FontSize',13)
ylabel('Frequncy','FontName','Times','FontSize',13)

print(gcf,'results/projection.png','-dpng','-r300')



%% APPENDIX: Numerical Flow Line model
function output = FlowLineM_init(geometry, parameter)
%%FLOWLINEM Flow line model, adapted from lab 7, Glaciology
%   Input
%       geometry:
%       parameter:
%
%   Output
%       output:

    % Geometry
    x = geometry.x;
    zb = geometry.zb;
    dx = geometry.dx;
    
    % Physical constants 
    rho_ice = 917; % kg/m^3, density of ice
    g       = 9.8; % m/s^2, gravitational constant
    n       = 3; % nonlinearity power
    A       = 3.8e-24; % Pa^-3 s^-1, flow law exponent
    secinyear = 365*24*3600; % seconds

    % Climate of the mountain
    % Different from before, dbdz is time-dependent (linear)
    dbdz = parameter.dbdz;   % Mass balance gradient (m/yr/m)
    dbdz = dbdz / secinyear; % Mass balance gradient (m/s/m)
    
    % ELA
    % similar to dbdz
    ELA = parameter.ELA;
    
    %
    % Define plot colors
    rockcolor = [0.7 0.5 0.2];   % in RGB space; it's light brown
    icecolor =  [0.7 0.9 1];     % in RGB space; it's light blue
    
    % initialize H
    H = 0;
    
    % Ice thickness in previous timestep -- to check if model has reached steady state
    H_previous = 0;   % meter
    L_previous = 0;

    % Initialize the flag that will stop the loop when model reaches steady state. 
    % Hint: does your model begin in steady state?
    steadystatereached = 0;

    % Initial surface slope: calculate using the diff function, 
    % and appending a slope of 0 at the beginning of the glacier. 
    dsdx = [0, diff(zb+H)/dx];   % m/m

    % Time counter: start at zero
    year = 0; % units: years

    % Plot the mountainside (no glacier on it yet)
    figure(1); clf; 
    subplot(2,1,1)
    area(x,zb,'facecolor',rockcolor,'basevalue',1000)

    % In the loop below, we will check for equilibrium every 10 years.
    last_year_checked = 0;  % Initialize the checkpoint of when we last made a plot and checked for equilibrium

    % Loop to run our model forward in time
    % Create a "while" loop that runs until the model has reached steady state. 
    % Use your variable called steadystatereached.
    while steadystatereached == 0 % when it is not at steady-state, keep crunching

        % Calculate the diffusivity, D. 
        % D has a different value at every x and D will change at every time step.
        D = 2*A*(rho_ice*g)^n.*H.^(n+2).*dsdx.^(n-1)/(n+1);   % m^2/s

        % Choose the biggest time step we can get away with (CFL condition)
        % but don't use a time step larger than 1 year
        dt = (dx^2)/(2*max(abs(D)));   % (This line is completed) unit is second
        dt = min(dt, secinyear);       

        % Calculate the ice flux AT each grid point
        F_pts  = 2*A.*(rho_ice*g*dsdx).^n.*H.^(n+2)/(n+1);   % m^2/s

        % We actually want to know the flux BETWEEN grid points.
        % We can achieve that by taking the mean of each pair of bordering grid
        % points. This is the flux at the midpoints.
        % (You won't have fluxes at the Nth point, because there is not a N+1th point to average the Nth point with.)
        % (So, F_mid will have N-1 elements, where N is the number of points in x.)
        F_left = F_pts(1:end-1);
        F_rigt = F_pts(2:end);
        F_mid = mean([F_left;F_rigt]);   % m^2/s

        % We weren't able to calculate the flux at the midpoints for the first and 
        % last point in the domain.  We want zero horizontal flux at these boundaries:
        % the only ice that goes into the first grid box is from snow; there is no ice 
        % flow into the top box. Same thing at the bottom box - there is no flow out 
        % the end of the glacier.
        % So we'll add boundaries with zero flux to get the flux, F, at each and every point in x.
        % (F will have N+1 elements, which includes the flux "between" the 0th and 1st point (0) and also the flux 
        % "between" the Nth and N+1th point (0).
        F = [0, F_mid, 0];   % m^2/s

        % We can now use F (flux) to calculate the rate of ice thickness change 
        % due to ice flowing in from upstream.  This is most of the RHS of the continuity equation, but
        % ignoring surface mass balance for now.
        dHdt = diff(F)/dx;   % m/s

        % We can also calculate the rate of ice thickness change from surface 
        % mass balance. You'll have to calculate the elevation of the glacier surface 
        % in order to apply the balance gradient. 
        Zs = H + zb; % m
        bdot = dbdz*(Zs-ELA);   % m/s
        %
        % Calculate the new ice thickness due to the two above changes that 
        % occurred over this time step:
        H = H + (bdot + dHdt)*dt;    % m

        % Set negative thicknesses to zero
        H(H < 0) = 0;  % where ice is gone

        % Calculate the new surface slope from the bed topography and the new ice thickness.
        % You should use the gradient function for this.
        dsdx = gradient(H+zb, dx);   % m/m

        % Advance the clock by the length of this time step (be mindful about units) 
        year = year + dt/secinyear;   % unit is year 

        % Check for steady state and re-plot the glacier geometry every time 10 years have passed
        if year - last_year_checked >= 10 % Check every 10 years

            figure(1);
            ax(1) = subplot(2,1,1,'replace'); hold on
            area(x/1e3,(zb+H)','facecolor',icecolor,'basevalue',5000);
            area(x/1e3,zb,'facecolor',rockcolor,'basevalue',5000)
            title(sprintf('Guoqu glacier with ELA at %1.0f m, year %1.0f',ELA,year))
            xlabel('Horizontal Distance (km)')
            ylabel('Elevation (m)')        
            % Now calculate u_surf at all grid points and plot it
            u = 2*A.*(rho_ice*g*abs(dsdx)).^n.*H.^(n+1)/(n+1);          % m
            ax(2) = subplot(2,1,2,'replace'); hold on
            plot(x/1e3,u*secinyear) % show as meters/year
            xlabel('Horizontal Distance (km)')
            ylabel('Velocity (m/yr)')
            %
            linkaxes(ax,'x')  % Allow the zoom tool to apply to both axes simultaneously
            %
            % Check whether the model has reached steady state. If the biggest change in ice thickness
            % anywhere on the glacier is less than 1 m, we'll consider the model to have reached steady state.
            L_now = x(find(H>0,1,'last'))/1e3;
            if max(H-H_previous) < 1 && L_now > 5 && abs(L_previous - L_now)< 0.5  % require 1 m precision to call it steady state
                fprintf('convergence criteria for steady state met after %1.0f years\n',year);
                % Update the variable steadystatereached so that your time loop won't have to run anymore
                steadystatereached = 1;

            elseif H(end) > 1 % if the glacier has grown beyond the domain
                % we leave the while loop
                break
            else% Model is not yet in steady state
                % Update the user on the time progress and glacier geometry
                fprintf('length %1.1fkm and max thickness change %1.1f m after %1.0f years\n',x(find(H>0,1,'last'))/1e3,max(H-H_previous),year)

                % Reset the year counter and model equilibrium variables
                H_previous = H;   
                L_previous = L_now; % (this line is complete)
                last_year_checked = floor(year);  % (this line is complete)
            end


        end

    end

    L = x(find(H>0,1,'last'));
    fprintf('Glacier length: %1.0f meters\n',L)
    
    % Output thickness
    output.H = H;
end
    
function output = FlowLineM_perturb(geometry, parameter, time)
%%FLOWLINEM Flow line model, adapted from lab 7, Glaciology
%   Input
%       geometry:
%       parameter:
%
%   Output
%       output:

    % Geometry
    x = geometry.x;
    zb = geometry.zb;
    dx = geometry.dx;
    %U_init = geometry.U_init;
    H = geometry.H_init;
    % change of climatic condition
    dbdz_0 = parameter.dbdz;   % Mass balance gradient (m/yr/m)
    ELA_0 = parameter.ELA;
    ddbdzdt = parameter.ddbdzdt;
    dELAdt = parameter.dELAdt;
    
    % Physical constants 
    rho_ice = 917; % kg/m^3, density of ice
    g       = 9.8; % m/s^2, gravitational constant
    n       = 3; % nonlinearity power
    A       = 3.8e-24; % Pa^-3 s^-1, flow law exponent
    secinyear = 365*24*3600; % seconds

    % Climate of the mountain
    % Different from before, dbdz is time-dependent (linear)
    t_temp = 0:0.01:time; % time axis for interp1
    dbdz_end = dbdz_0 + ddbdzdt*time;
    dbdz_t = interp1([0,time], [dbdz_0, dbdz_end], t_temp);
    dbdz_t = dbdz_t / secinyear; % Mass balance gradient (m/s/m)
    
    % ELA
    % similar to dbdz
    ELA_end = ELA_0 + dELAdt*time;
    ELA_t = interp1([0,time], [ELA_0, ELA_end], t_temp);

    %
    % Define plot colors
    rockcolor = [0.7 0.5 0.2];   % in RGB space; it's light brown
    icecolor =  [0.7 0.9 1];     % in RGB space; it's light blue
    %
    % Initial surface slope: calculate using the diff function, 
    % and appending a slope of 0 at the beginning of the glacier. 
    dsdx = [0, diff(zb+H)/dx];   % m/m
    
    % previous H
    H_previous = 0;

    % Time counter: start at zero
    year = 0; % units: years

    % Plot the mountainside (no glacier on it yet)
    figure(1); clf; 
    subplot(2,1,1)
    area(x,zb,'facecolor',rockcolor,'basevalue',1000)

    % In the loop below, we will check for equilibrium every 10 years.
    last_year_checked = 0;  % Initialize the checkpoint of when we last made a plot and checked for equilibrium
    
    % time-integrated accumulation over the year
    accum = 0; % meter

    % Loop to run our model forward in time
    % Create a "while" loop that runs until the model has reached steady state. 
    % Use your variable called steadystatereached.
    while year < time-1 % when it is not at steady-state, keep crunching

        % Calculate the diffusivity, D. 
        % D has a different value at every x and D will change at every time step.
        D = 2*A*(rho_ice*g)^n.*H.^(n+2).*dsdx.^(n-1)/(n+1);   % m^2/s

        % Choose the biggest time step we can get away with (CFL condition)
        % but don't use a time step larger than 1 year
        dt = (dx^2)/(2*max(abs(D)));   % (This line is completed) unit is second
        dt = min(dt, secinyear);       

        % Calculate the ice flux AT each grid point
        F_pts  = 2*A.*(rho_ice*g*dsdx).^n.*H.^(n+2)/(n+1);   % m^2/s

        % We actually want to know the flux BETWEEN grid points.
        % We can achieve that by taking the mean of each pair of bordering grid
        % points. This is the flux at the midpoints.
        % (You won't have fluxes at the Nth point, because there is not a N+1th point to average the Nth point with.)
        % (So, F_mid will have N-1 elements, where N is the number of points in x.)
        F_left = F_pts(1:end-1);
        F_rigt = F_pts(2:end);
        F_mid = mean([F_left;F_rigt]);   % m^2/s

        % We weren't able to calculate the flux at the midpoints for the first and 
        % last point in the domain.  We want zero horizontal flux at these boundaries:
        % the only ice that goes into the first grid box is from snow; there is no ice 
        % flow into the top box. Same thing at the bottom box - there is no flow out 
        % the end of the glacier.
        % So we'll add boundaries with zero flux to get the flux, F, at each and every point in x.
        % (F will have N+1 elements, which includes the flux "between" the 0th and 1st point (0) and also the flux 
        % "between" the Nth and N+1th point (0).
        F = [0, F_mid, 0];   % m^2/s

        % We can now use F (flux) to calculate the rate of ice thickness change 
        % due to ice flowing in from upstream.  This is most of the RHS of the continuity equation, but
        % ignoring surface mass balance for now.
        dHdt = diff(F)/dx;   % m/s

        % We can also calculate the rate of ice thickness change from surface 
        % mass balance. You'll have to calculate the elevation of the glacier surface 
        % in order to apply the balance gradient. 
        Zs = H + zb; % m
        % interpolate to get the new ELA and dbdz value
        ELA_now  = interp1(t_temp, ELA_t, year + dt/secinyear);
        dbdz_now = interp1(t_temp, dbdz_t, year + dt/secinyear);
        bdot = dbdz_now*(Zs-ELA_now);   % m/s
        %
        % Calculate the new ice thickness due to the two above changes that 
        % occurred over this time step:
        H = H + (bdot + dHdt)*dt;    % m

        % Set negative thicknesses to zero
        H(H < 0) = 0;  % where ice is gone

        % Calculate the new surface slope from the bed topography and the new ice thickness.
        % You should use the gradient function for this.
        dsdx = gradient(H+zb, dx);   % m/m

        % Advance the clock by the length of this time step (be mindful about units) 
        year = year + dt/secinyear;   % unit is year 
        
        % calculate integrated accumulation over the years
        b_core = interp1(Zs, bdot*dt, 5750); 
        accum = accum + b_core;

    end

    L = x(find(H>0,1,'last'));
    fprintf('Glacier length: %1.0f meters\n',L)
    
    % total melt volume
    melt = sum(H - geometry.H_init);
    
    % Output thickness
    output.H = H;
    output.L = L;
    output.accum = accum;
    output.melt = melt;
end
    
function output = FlowLineM_proj(geometry, parameter, time)
%%FLOWLINEM Flow line model, adapted from lab 7, Glaciology

    % Geometry
    x = geometry.x;
    zb = geometry.zb;
    dx = geometry.dx;
    %U_init = geometry.U_init;
    H = geometry.H_init;
    % change of climatic condition
    dbdz_0 = parameter.dbdz;   % Mass balance gradient (m/yr/m)
    ELA_0 = parameter.ELA;
    ddbdzdt = parameter.ddbdzdt;
    dELAdt = parameter.dELAdt;
    
    % Physical constants 
    rho_ice = 917; % kg/m^3, density of ice
    g       = 9.8; % m/s^2, gravitational constant
    n       = 3; % nonlinearity power
    A       = 3.8e-24; % Pa^-3 s^-1, flow law exponent
    secinyear = 365*24*3600; % seconds

    % Climate of the mountain
    % Different from before, dbdz is time-dependent (linear)
    t_temp = 0:0.01:time; % time axis for interp1
    dbdz_end = dbdz_0 + ddbdzdt*time;
    dbdz_t = interp1([0,time], [dbdz_0, dbdz_end], t_temp);
    dbdz_t = dbdz_t / secinyear; % Mass balance gradient (m/s/m)
    
    % ELA
    % similar to dbdz
    ELA_end = ELA_0 + dELAdt*time;
    ELA_t = interp1([0,time], [ELA_0, ELA_end], t_temp);

    % Initial surface slope: calculate using the diff function, 
    % and appending a slope of 0 at the beginning of the glacier. 
    dsdx = [0, diff(zb+H)/dx];   % m/m
    
    % previous H
    H_previous = 0;

    % Time counter: start at zero
    year = 0; % units: years
    
    % time-integrated accumulation over the year
    accum = 0; % meter
    
    % create time axis for recording
    time_ax = [];
    dmelt = [];

    % Loop to run our model forward in time
    % Create a "while" loop that runs until the model has reached steady state. 
    % Use your variable called steadystatereached.
    while year < time-1 % when it is not at steady-state, keep crunching

        % Calculate the diffusivity, D. 
        % D has a different value at every x and D will change at every time step.
        D = 2*A*(rho_ice*g)^n.*H.^(n+2).*dsdx.^(n-1)/(n+1);   % m^2/s

        % Choose the biggest time step we can get away with (CFL condition)
        % but don't use a time step larger than 1 year
        dt = (dx^2)/(2*max(abs(D)));   % (This line is completed) unit is second
        dt = min(dt, secinyear);       

        % Calculate the ice flux AT each grid point
        F_pts  = 2*A.*(rho_ice*g*dsdx).^n.*H.^(n+2)/(n+1);   % m^2/s

        % We actually want to know the flux BETWEEN grid points.
        % We can achieve that by taking the mean of each pair of bordering grid
        % points. This is the flux at the midpoints.
        % (You won't have fluxes at the Nth point, because there is not a N+1th point to average the Nth point with.)
        % (So, F_mid will have N-1 elements, where N is the number of points in x.)
        F_left = F_pts(1:end-1);
        F_rigt = F_pts(2:end);
        F_mid = mean([F_left;F_rigt]);   % m^2/s

        % We weren't able to calculate the flux at the midpoints for the first and 
        % last point in the domain.  We want zero horizontal flux at these boundaries:
        % the only ice that goes into the first grid box is from snow; there is no ice 
        % flow into the top box. Same thing at the bottom box - there is no flow out 
        % the end of the glacier.
        % So we'll add boundaries with zero flux to get the flux, F, at each and every point in x.
        % (F will have N+1 elements, which includes the flux "between" the 0th and 1st point (0) and also the flux 
        % "between" the Nth and N+1th point (0).
        F = [0, F_mid, 0];   % m^2/s

        % We can now use F (flux) to calculate the rate of ice thickness change 
        % due to ice flowing in from upstream.  This is most of the RHS of the continuity equation, but
        % ignoring surface mass balance for now.
        dHdt = diff(F)/dx;   % m/s

        % We can also calculate the rate of ice thickness change from surface 
        % mass balance. You'll have to calculate the elevation of the glacier surface 
        % in order to apply the balance gradient. 
        Zs = H + zb; % m
        % interpolate to get the new ELA and dbdz value
        ELA_now  = interp1(t_temp, ELA_t, year + dt/secinyear);
        dbdz_now = interp1(t_temp, dbdz_t, year + dt/secinyear);
        bdot = dbdz_now*(Zs-ELA_now);   % m/s
        %
        % Calculate the new ice thickness due to the two above changes that 
        % occurred over this time step:
        H = H + (bdot + dHdt)*dt;    % m

        % Set negative thicknesses to zero
        H(H < 0) = 0;  % where ice is gone

        % Calculate the new surface slope from the bed topography and the new ice thickness.
        % You should use the gradient function for this.
        dsdx = gradient(H+zb, dx);   % m/m

        % Advance the clock by the length of this time step (be mindful about units) 
        year = year + dt/secinyear;   % unit is year 
        
        % append time axis and volume change
        time_ax = [time_ax, year];
        if H_previous == 0 % first iteration; we just let dmelt be 0
            dmelt = [dmelt, 0];
        else
            dmelt = [dmelt, sum(H - H_previous)*dx];
        end
        
        % update H_previous
        H_previous = H;


    end

    L = x(find(H>0,1,'last'));
    fprintf('Glacier length: %1.0f meters\n',L)
    
    % total melt volume
    tot_melt = sum(H - geometry.H_init);
    
    % Output thickness
    output.H = H;
    output.L = L;
    output.accum = accum;
    output.tot_melt = tot_melt;
    output.time_ax = time_ax;
    output.dmelt = dmelt;
end
