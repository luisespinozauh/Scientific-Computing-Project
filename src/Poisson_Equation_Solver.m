% Scientific Computing
% Luis Espinoza, 1226327
% Project A - Poisson Equation
% Poisson Equation Solver using the Gauss- Seidel Method

clc
clear all

%% Define dimension of 2-D grid

a_x=0;
a_y=0;
b_x= (2*pi);  
b_y= (2*pi);

N=10;    % Number of Grid Points, change for grid convergence study  
Nx=N;    % Number of nodes in x-direction
Ny=N;    % Number of nodes in y-direction
tol=1e-06;   % Tolerance
err= 1;      % Error
iter=0;      % Iteration counter
frequency=10; % how often to perform checkpointing (intermediate save)

x=linspace(a_x, b_x, Nx);   % Mesh
y=linspace(a_y, b_y, Ny);
h=x(2)-x(1);               % Step Size

%% Define u

u=zeros(Nx,Ny);

u(Nx,:)=x.*((2*pi)-x).^2;             % Given Boundary Conditions, optimization- instead of using for loops
u(1,:)=(((2*pi)-x).^2).*cos(x/2);
u(:,1)=(4*pi*pi)-((2*pi).*y);
u(:,Ny)=b_x;

F=sin(x/(2*pi))'*cos((y+pi)/2);      % Forcing function, optimization- instead of using a for loop

%% Checkpoint

% Before the start of the iteration loop, "check-in" each variable
% that should be checkpointed in the event of restarting the job

matfile = 'PoissonEquationSolution.mat';     % mandatory; name of checkpoint mat-file
s = struct();               % mandatory; create struct for checkpointing
s = chkin(s,{'iter'});      % mandatory; iter is iteration loop index
s = chkin(s,{'frequency'}); % mandatory; frequency is checkpointing period 
                            % i.e., how often to perform a save

% continue until all variables are checked in. Note that you are only
% checking in the variables, they don't need to have been already defined

chkNames = fieldnames(s);    % the full list of variables to checkpoint
nNames = length(chkNames);   % number of variables in list

%% Iterative Looping- Gauss-Seidel Method

tic;        % Timer to evaluate Performance

while max(err(:)) > tol
    iter= iter + 1;
    %pause(1);       % To test the restart script  
    
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

timedoc=toc;
fprintf('Number of iterations is %f.\n',iter)
fprintf('Running Time is %f seconds.\n',timedoc)

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