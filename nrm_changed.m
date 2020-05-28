%
%  nrm_changed.m	Run Next Reaction Method on the model
%
%                   X1 + X2 -> X3, X3 -> 0 
%
%               with rates c1 = 0.001, c2 = 0.001.  X1 + X2 -> X3 
%               affects X1 and X2 at initiation and X3 at completion.
%
%  usage:	nrm_changed(y,maxi,tau,tend,runs)
%
%		y = initial # of [A B] molecules
%
%		maxi = number of iterations to be taken
%		
%               tau = delay
%
%               tend = max time for the program to run
%
%               runs = number of times to run to get stats.
%
%  example:	nrm_changed([1000 1000 0],10000, 0.1, 1, 10)
%

function [Time CPUtime] = nrm_changed(y,maxi,tau,tend,runs)

  
format compact
rand('state',sum(100*clock))

bt = clock;

% count will let me know hom many iterations took place before my break
%count = 0;

Delta = 0;

%These generate clock times.

%clocktime = zeros([runs,1]);

%rates.
c = [0.001 0.001];

%storage for stats
%L = [runs,3];

for j = 1:runs

%ct = clock;

%data storage for delayed reaction times

s = [inf];
Td = 1;

%initialize.
T = 0;%zeros([maxi,1]);
x = zeros([maxi,3]);
a = zeros([2,1]);

%initialize the x values
x(1,:) = y;

%Set Poisson times and next Poisson Jumps.  Also set t_ks.
Tk = zeros([2,1]);
Pk = zeros([2,1]);
t = zeros([2,1]);

r1 = rand;
r2 = rand;

Pk(1) = log(1/r1);
Pk(2) = log(1/r2);

for i=1:maxi-1,
 
  if T >= tend
    break
  end
  if sum(x(i,:)) == 0
    break
  end

  a(1) = c(1)*x(i,1)*x(i,2);
  a(2) = c(2)*x(i,3);

  if a(1) ~= 0
    t(1) = (Pk(1) - Tk(1))/a(1);
  else
    t(1) = inf;
  end

  if a(2) ~= 0
    t(2) = (Pk(2) - Tk(2))/a(2);
  else
    t(2) = inf;
  end 

  t0 = s(1);

  if t0 < t(1) && t0 < t(2)          %apply completion of delayed reaction

    Delta = t0;
    T = T + Delta;
    x(i+1,1) = x(i,1);
    x(i+1,2) = x(i,2);  
    x(i+1,3) = x(i,3)-1;

    s = s(1,2:Td);
    Td = Td - 1;
    s = s - Delta;

  elseif t(1) < t(2) %apply reaction 1

    Delta = t(1);
    T = T + Delta;
    x(i+1,1) = x(i,1) - 1;
    x(i+1,2) = x(i,2) - 1;  
    x(i+1,3) = x(i,3) + 1;

    s = s - Delta;

    r1 = rand;
    Pk(1) = Pk(1) + log(1/r1);

  else                        %apply reaction 2 and store

    Delta = t(2);
    T = T + Delta;
    x(i+1,1) = x(i,1);
    x(i+1,2) = x(i,2);  
    x(i+1,3) = x(i,3);

    s = s - Delta;

    s = [s(1:Td - 1), tau, inf];    
    Td = Td + 1;
    r2 = rand;
    Pk(2) = Pk(2) + log(1/r2);

  end % if t0 < t(1) && t0 < t(2) %apply reaction 1

  Tk(1) = Tk(1) + a(1)*Delta;
  Tk(2) = Tk(2) + a(2)*Delta;


end % for i=1:maxi-1,

%for b = 1:3
%L(j,b) = x(i,b);
%end

%cputime(j) = etime(clock,ct);

end %for j = 1:runs

CPUTime = etime(clock,bt)

%     count = i;

%     X1 = mean(L(:,1));
%     X2 = mean(L(:,2));
%     X3 = mean(L(:,3));

%     Time = T(count);


%The following is so I can print out the clock times.
%Index = zeros([runs,1]);

%for b = 1:runs-1
%  Index(b+1) = Index(b) + 1;
%end

%figure
%plot(Index,cputime,'b')
%xlabel('Simulation Step')
%ylabel('CPU time for timestep')

%figure
%plot(T(1:count),x(1:count,1),'b',T(1:count),x(1:count,2),'r',T(1:count),x(1:count,3),'g')
%grid on
%legend('X_1','X_2','X_3')
%xlabel('Simulation Time t')
%ylabel('# of Molecules')


end %The program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
