% Luis Espinoza, 1226327
% Restart Script for SOR Method after a system crash
% Sometimes files take a long time to run to completion. As a result, sometimes they crash due to a variety of reasons: power failure, walltime limit, scheduled shutdown, etc. 
% Checkpoint/Restarting has long been a common technique to tackle this issue. Checkpointing/Restarting essentially means saving data to disk periodically so that, if need be, 
% you can restart the job from the point at which your data was last saved.                          %


%% Any information defined in original code that has not been checkpointed should be placed here.


%% Restarting
matfile = 'PoissonEquationSolution';   % should match that in test_checkpoint.m
load(matfile);        % retrieve data from matfile
iter1 = iter+1;       % iter is the last time test_checkpoint issued
                      % a save; we start computing on the next step

%% Iterative Looping- SOR                  

while max(err(:)) > tol
    iter= iter + 1;
    
    % If you want to test the restart script, use the function pause(1) to slow down the while loop.
    % This will slow down the while loop to 1 sec per iteration so that ctrl + C can used be to
    % "kill" the code to simulate a computer crash. From there, use the restart script to restart the loop.  
    
    if mod(iter, frequency) == 0 % If statement, checkpoints periodically (determined by the frequency)
        chkpt                    % chkpt script performs checkpointing (save) every *frequency* iterations
        fprintf(1, ['Checkpointing frequency is every %2d iterations.' ...
          'Data updated at iteration %3d\n'], ...
          frequency, iter);      % Confirm after each checkpointing event 
    end
    
    uold=u;
    
    for i=2:Nx-1        % Gauss-Seidel Method will only be performed on interior nodes
        for j=2:Ny-1
            u(i,j)=(B.*(1/4)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)+(h.^2)*F(i,j)))+(1-B)*u(i,j); % Discretizatized Poisson Equation
        end
        u(i,Nx)=(B.*(1/4)*(u(i+1,j)+u(i+1,j)+u(i,j+1)+u(i,j-1)+(h.^2)*F(i,j)))+(1-B)*u(i,j);    % Neumann Boundary Condition
    end
    u(N,N)=0;
    u(1,N)=0;
    unew=u;
    err=abs((uold-unew)./unew);     % Equation for relative error. We only consider the maximum error of all interior nodes
    fprintf(1, 'Completed iteration %d\n', iter);
end

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

