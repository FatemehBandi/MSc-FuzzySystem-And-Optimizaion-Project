
%___________________________________________________________________%
%                                                                   %
%  Developed in MATLAB R2014a                                       %
%     							                                          	    %
%								                                                    %
%                                                                   %
%  programmer: Fatemeh Bandi		                                    %
%								                                                    %
%								                                                    %
%        paper: 	                              				            %
%        Novel adaptive hybrid rule network based on TS fuzzy rules %
%        using an improved quantum-behaved particle swarm 	        %
%        optimization for medical data classification               %
%                                                                   %
%								                                                    %
%___________________________________________________________________%




clc;
clear ;
close all;

%% load data
load mgdata.dat
x  =  mgdata(:, 2);

for t=118:1117
inputData(t-117,:) = [x(t-18) x(t-12) x(t-6) x(t) ];
outputData(t-117,:) = [x(t+6)];
end

% train_data & test_data
trninput = inputData(1:500,:);
chkinput = inputData(501:end,:);

trnoutput = outputData(1:500,:);
chkoutput = outputData(501:end,:);

ntrnData = size(trninput,1);
nchkData = size(chkinput,1);

%% calculate  Intervals
  global CI
  Ri = [2 2 2 2];
  CI = [Ri(1)+1 Ri(1)+1 Ri(1)+1 Ri(1)+1];
   
  nci = numel(Ri);
  RI = zeros(nci,max(Ri)+1);
  maxt = max(max(trninput));
  mint = min(min(trninput));
  A = zeros(nci,max(CI)+1);
  

  % calculate RI
  
  for i=1:nci
      
      Ris(i) = (maxt-mint)/Ri(i);
      RI(i,1) = mint;
      for j=1:Ri(i)
          RI(i,j+1) = RI(i,j)+Ris(i); % CIs is size of each interval 
      end
      
  end

  % calculate CI
  
  for i=1:nci
      
      A(i,1) = RI(i,1)-((RI(i,1)+RI(i,2))/2);
      for j=2:CI(i)
        
          A(i,j) = (RI(i,j-1)+RI(i,j))/2;
      end
      A(i,j+1) = RI(i,j)+((RI(i,j-1)+RI(i,j))/2);
  end
  
  
  %% mean &spread 
  mo = zeros(nci,max(CI));
  spread = nan(nci,1);
  for i=1:nci
     for j=1:CI(i)
     mo(i,j) = (A(i,j)+A(i,j+1))/2;
     
     end
     spread(i) = mo(i,1)-A(i,1);
  end
  
  %% initialization
  
  npop = 30; %number of particle
  MAXITER = 10;
  Err = 0.01;
  beta = .3;
  alpha = .4;
  
  
  % lb,ub [w11,...,wnn,a11,...,ann,b1,...,bn,RL,Rs]
  lb = zeros(1,(numel(CI)*CI(1))*2+CI(1)+2);
  ub = zeros(1,(numel(CI)*CI(1))*2+CI(1)+2);
 
  
  
 %%%%%%%% calculate rang W
 
 % first phase: calculate how many data is belong to each CI
  
 v = zeros(nci,max(CI));
 W = zeros(nci,max(CI));
 for i=1:nci
     for j=1:CI(i)
         for k=1:500
          if A(i,j) <= trninput(k,i) & trninput(k,i)<A(i,j+1)
              v(i,j) = v(i,j)+1;
          end
         end
     end
     
 end
  
 % second phase:
 
  for i=1:nci
     for j=1:CI(i)
       s=sum(v,2)';
       W(i,j) = v(i,j)/s(i);
     end
  end
  
  
%%%% Range of lower & upper bound  
      cn=0;
   
  for j=1:numel(CI)
      n=0;
      for i=1:CI(1)-1
          lb(i+1+((j-1)*cn)) = W(j+n);
          n = numel(CI);
      end
      n1=0;
      for i=1:CI(1)
          ub(i+((j-1)*cn)) = W(j+n1);
          n1 = numel(CI);
      end
      
      cn=CI(1);
  end
  
  
  t=numel(CI)*CI(1);
  for l=t+1:2*t
     lb(l) = 0;
     ub(l) = 1;
  end
 
  for k=l+1:size(lb,2)-2
     lb(k) = -3;
     ub(k) = 1;
  end
  
  lb(k+1) = 1;
  ub(k+1) = CI(1);
  lb(k+2) = 0;
  ub(k+2) = 1;
  
  % creat population
  
  ind.position=[];
  ind.cost=[];
  dimension=size(lb,1);
  particle=repmat(ind,npop,1);
  nvar=size(lb,2); %nvar is number of gen
  
  for i=1:npop
     particle(i).position=unifrnd(lb,ub,1,nvar);
     y=membersh( particle(i).position,A,mo,spread,trninput);
     y=y';
     particle(i).cost=COST(y,trnoutput);
  end

  pParticle=particle; %pbest
  gbest=zeros(1,dimension);
  
  [v, index]=min([particle.cost]);
  gbest=particle(index);%gbest
  
  %% main QPSO
  
  for t=1:MAXITER
  
       m = 0;
      start_time_train=cputime;% start CPU time
      
      for o=1:npop
       m = m + (minus(pParticle(o).position,gbest.position));
       
      end
    
      mbest=alpha.*gbest.position+((1-alpha)/(npop-1)).*m; 
      
      for i=1:npop  
        k=rand();
        u=rand(1,dimension);
        b=beta.*(mbest-particle(i).position);
        v=-log(u);
        if k>=0.5
        yhat=pParticle(i).position+b.*v;
         else
         yhat=pParticle(i).position-b.*v;   
        end
         particle(i).position=yhat;
       
          
         y=membersh( particle(i).position,A,mo,spread,trninput);
         y=y';
         particle(i).cost=COST( y,trnoutput);
         
            if particle(i).cost<pParticle(i).cost
                pParticle(i).position=particle(i).position;
                pParticle(i).cost=particle(i).cost;
            end
            if pParticle(i).cost<gbest.cost
                gbest.position=pParticle(i).position;
                gbest.cost=pParticle(i).cost;
            end
      end
      
      
     y=membersh(gbest.position,A,mo,spread,trninput);
     y=y';
     Error(t)=COST(y,trnoutput);
     
     if Error(t)<Err
         break;
     end
     
     end_time_train=cputime;
     TrainingTime(t)=end_time_train-start_time_train ;       %   Calculate CPU time (seconds) spent for training iteration

     disp(['it('  num2str(t) ')'   'Error='  num2str( Error(t)) ]);
     
     if Error(t)<Err
         break;
     end
     
  end
  
  
  avg_Time=sum(TrainingTime)/t
  
  % plot Error & Iteration
  semilogy(Error,'b:');
  xlabel('iteration');
  ylabel('Cost');
  
  % test Error
  
   y=membersh(gbest.position,A,mo,spread,chkinput);
   y=y';
   ErrorTest=COST(y,chkoutput);
   disp([  'Error Test='  num2str( ErrorTest) ]);
