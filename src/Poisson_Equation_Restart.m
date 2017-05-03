%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On a cluster such as SCC, jobs that take hours or days to run to completion    %
% on occasion could suffer from system abort in mid-stream due to a variety of   %
% reasons: power failure, walltime limit, scheduled shutdown, to name a few.     %
% Rerunning jobs oftens means inevitable delay and waste of system resources.    %
% Checkpoint Restarting has long been a common technique with which researchers  %
% turn to tackle this issue. In brief, Checkpoint Restarting essentially means   %
% saving data to disk periodically so that if need be, you can restart the job   %
% where data were last saved. Checkpoint Restarting can either be dealt with     %
% through a batch scheduler (if supported) or do it the old-fashioned way by     %
% writing to disk manually yourself. On the SCC, the OGE batch scheduler         %
% checkpointing feature is not supported. A rather simple MATLAB checkpointing   %
% tool has been developed by SCV to ease the pain of manual checkpointing.       %
% This MATLAB example code demonstrates the usage of this checkpointing tool to  %
% facilitate rerunning of your job should a system abort occurs.                 %
% The problem is to compute the arithmetic                                       %
% A = 1 + 2 + 3 + . . . + N = N*(N+1)/2                                          %
%                                                                                %
% December 8, 2013                                                               %
% Kadin Tseng, SCV, Boston University  (kadin@bu.edu)                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Any information defined in original code that has not
%%% been checkpointed should be placed here.

%===============================================================
matfile = 'PoissonEquationSolution';   % should match that in test_checkpoint.m
load(matfile);        % retrieve data from matfile
iter1 = iter+1;       % iter is the last time test_checkpoint issued
                      % a save; we start computing on the next step
%===============================================================
while max(err(:)) > tol
    iter= iter + 1;
     
            
    % Checkpoints periodically (determined by the determined frequency)
    if mod(iter, frequency) == 0
        chkpt   % performs checkpointing (save) every *frequency* iterations
        fprintf(1, ['Checkpointing frequency is every %2d iterations.' ...
          'Data updated at iteration %3d\n'], ...
          frequency, iter);  % Confirm after each checkpointing event
      
    end
    
    uold=u;
    
    for i=2:Nx-1
        for j=2:Ny-1
            u(i,j)=0.25*(u(i+1,j)+ u(i-1,j)+ u(i,j+1)+ u(i,j-1)+ (F(i,j)*(h^2)));
        
        end
    end
    
    unew=u;
    err=abs((uold-unew)./unew);     % Equation for relative error. We only consider the maximum error
    fprintf(1, 'Completed iteration %d\n', iter);
end

timedoc2=toc;
fprintf('Number of iterations is %f.\n',iter)
fprintf('Running Time is %f seconds.\n',timedoc2)

%% Plot the results

figure
contourf(u)
h=gca; 
get(h,'FontSize') ;
set(h,'FontSize',12);
colorbar('location','eastoutside','fontsize',12);
axis([1 Nx 1 Ny]);
xlabel('X (Number of Nodes in X-direction)','fontSize',12);
ylabel('Y (Number of Nodes in Y-direction)','fontSize',12);
title('Numerical Solution to Poisson Equation, u(x,y)','fontsize',12);
fh = figure(1);
set(fh, 'color', 'white');

figure
mesh(x,y,u)
xlabel('X domain','fontSize',12);
ylabel('Y domain','fontSize',12);
zlabel('u(x,y)')
title('Numerical Solution to Poisson Equation','fontsize',12);

