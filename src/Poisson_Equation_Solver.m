% Scientific Computing Class
% Luis Espinoza, 1226327
% Project A - Poisson Equation
% Poisson Equation Solver using the Gauss- Seidel Method

clc
clear all
close all
%% Define dimension of 2-D grid

a_x=0;   % Domain and Range
a_y=0;
b_x= (2*pi);  
b_y= (2*pi);

N=100;    % Number of Grid Points, change for grid convergence study  
Nx=N;    % Number of nodes in x-direction
Ny=N;    % Number of nodes in y-direction
tol=1e-06;    % Tolerance
err= 1;       % Error
iter=0;       % Iteration counter
frequency=10; % how often to perform checkpointing (intermediate save)

x=linspace(a_x, b_x, Nx);   % Mesh
y=linspace(a_y, b_y, Ny);
h=x(2)-x(1);                % Step Size

%% Define u

u=zeros(Nx,Ny);                       % Initial guess for Gauss-Seidel 
                                      % Method is zero for all interior nodes

u(Nx,:)=x.*((2*pi)-x).^2;             % Given Boundary Conditions (3 Dirichlet and 1 Neumann)
u(1,:)=(((2*pi)-x).^2).*cos(x/2);     % Optimized code by using vectorization instead of for loop
u(:,1)=(4*pi*pi)-((2*pi).*y);           
                            

F=sin(x/(2*pi))'*cos((y+pi)/2);       % Forcing function, optimization- instead of using a for loop

%% Checkpointing
% Sometimes files take a long time to run to completion. As a result, sometimes they crash due to a variety of reasons: power failure, walltime limit, scheduled shutdown, etc. 
% Checkpoint/Restarting has long been a common technique to tackle this issue. Checkpointing/Restarting essentially means saving data to disk periodically so that, if need be, 
% you can restart the job from the point at which your data was last saved. 


% Before the start of the iteration loop, "check-in" each variable
% that should be checkpointed in the event of restarting the job

matfile = 'PoissonEquationSolution.mat';     % mandatory; name of checkpoint mat-file
s = struct();                                % mandatory; create struct for checkpointing
s = chkin(s,{'iter'});                       % mandatory; iter is iteration loop index
s = chkin(s,{'frequency'});                  % mandatory; frequency is checkpointing period 
                                             % i.e., how often to perform a save

% continue until all variables are checked in. Note that you are only
% checking in the variables, they don't need to have been already defined

chkNames = fieldnames(s);    % the full list of variables to checkpoint
nNames = length(chkNames);   % number of variables in list

%% Iterative Looping- Gauss-Seidel Method

tic;        % Timer to evaluate Performance

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
            u(i,j)=0.25*(u(i+1,j)+ u(i-1,j)+ u(i,j+1)+ u(i,j-1)+ (F(i,j)*(h^2))); % Discretizatized Poisson Equation
        end
        u(i,Nx)=0.25*(u(i+1,j)+ u(i+1,j)+ u(i,j+1)+ u(i,j-1)+ (F(i,j)*(h^2)));     % Neumann Boundary Condition
    end
    
    unew=u;
    err=abs((uold-unew)./unew);     % Equation for relative error. We only consider the maximum error of all interior nodes
    fprintf(1, 'Completed iteration %d\n', iter);
end

timedoc=toc;
fprintf('Number of iterations is %f.\n',iter)
fprintf('Running Time is %f seconds.\n',timedoc)

%% Plot the results

figure
contourf(x,y,u)                    % Contour Plot
h=gca; 
get(h,'FontSize') ;
set(h,'FontSize',12);
colorbar('location','eastoutside','fontsize',12);
xlabel('X (Number of Nodes in X-direction)','fontSize',12);
ylabel('Y (Number of Nodes in Y-direction)','fontSize',12);
title('Numerical Solution to Poisson Equation, u(x,y)','fontsize',12);
fh = figure(1);
set(fh, 'color', 'white');

figure
mesh(x,y,u)                    % Surface Plot
xlabel('X domain','fontSize',12);
ylabel('Y domain','fontSize',12);
zlabel('u(x,y)')
title('Numerical Solution to Poisson Equation','fontsize',12);
