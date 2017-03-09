function [epsMax, epsdotMax] = strainFromDisplacement(fname, numTimeSteps)
% Strain from AFiNeS simulation outputs measured as average displacement
% normalized to the linear dimension of the simulation box, and strain rate 
% measured as a simple forward difference of the strain
% 
% INPUTS    fname        : String, name of file to be analyzed
%           numTimeSteps : Scalar, number of time steps to measure strain 
%                          (i.e. numTimeSteps=2 measures positions on every other frame available)
% 
% OUTPUTS   Saves the following data and plots in fname folder:
% 
%           eps          : Time series of strain of the network, averaged at each time step over all the beads and
%                          normalized to the linear dimension of the simulation box, specifically:
%                              eps(t) = <|x(t)-x(0)|>/L
%                          where L is calculated using the position data in simdata.mat
%           epsdot       : Time series of the strain rate of the network, specifically:
%                              epsdot(t) = (eps(t+dt)-eps(t))/dt   
%               
%           Returns the following variables:
%               
%           epsMax       : Maximum of eps
%           epsdotMax    : Maximum of epsdot
%        
% Assumes:  -Periodic Boundary Conditions
%           -There exists a data file called "simdata.mat" that contains the position vs time data of the actin beads
%
% Created by Daniel Seara at 2017/01/30 13:20
    tic;

    cd(fname);
    disp(fname)

    run_read_data(fname)

    dataFile = 'simdata.mat';
    vars = {'adata', 'params'};
    load(dataFile, vars{:});

    % Find time step between frames
    DeltaT = numTimeSteps*params.dt;

    % Enforce periodic boundary conditions
    xyDisplacement = mod(diff(adata(:,1:2, 1:numTimeSteps:end),1,3) + params.L/2, params.L) - params.L/2;
    dx = squeeze(xyDisplacement(:,1,:)); cumDX = cumsum(dx,2);
    dy = squeeze(xyDisplacement(:,2,:)); cumDY = cumsum(dy,2);
    cumDR = [zeros(size(dx,1),1), sqrt(cumDX.^2 + cumDY.^2)]./params.L;

    eps = mean(cumDR);
    epsMax = max(eps);

    plot(params.timestep,eps);
    xlabel('time (s)')
    ylabel('$\epsilon (t)$', 'Interpreter', 'latex')
    title(fname, 'Interpreter', 'none')
    
    % Make figure prettier...
    ax=findall(gca,'Type','line');
    for i=1:length(ax)
        set(ax(i),'linewidth',2);
    end
    ax=findall(gcf,'Type','text');
    for i=1:length(ax)
        set(ax(i),'fontname','times','fontsize',14);
    end
    saveas(gcf, 'strainSeries', 'fig');
    saveas(gcf, 'strainSeries', 'epsc');
    
    clf reset;

    % epsdot = diff(eps)./params.dt;
    % epsdotMax = max(epsdot);

    % plot(timestep(2:end-1),epsdot(2:end)); %first data point of epsdot seems to be huge all the time? Unsure why...
    % xlabel('time (s)')
    % ylabel('$\dot{\epsilon} (t) \; (s^{-1})$', 'Interpreter', 'latex')
    % title(fname, 'Interpreter', 'none')
    
    % % Make figure prettier...
    % ax=findall(gca,'Type','line');
    % for i=1:length(ax)
    %     set(ax(i),'linewidth',2);
    % end
    % ax=findall(gcf,'Type','text');
    % for i=1:length(ax)
    %     set(ax(i),'fontname','times','fontsize',14);
    % end
    % saveas(gcf, 'strainRateSeries', 'fig');
    % saveas(gcf, 'strainRateSeries', 'epsc');

    % clf reset;

    % clearvars -except eps epsdot params epsMax epsdotMax

    % if ispc
    %     save([pwd,'\straindata.mat'], 'eps', 'epsdot', 'params');
    % elseif isunix
    %     save([pwd,'/straindata.mat'], 'eps', 'epsdot', 'params');
    % end

    cd ..

    % clearvars -except epsMax epsdotMax
    toc
end % end of function