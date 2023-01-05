function y=rules(x,ind1,WF,a,RS,RL,SI)
global CI

     if RL==1
       SI1=zeros(numel(CI),RL);
       for i=1:numel(CI)
           SI1(i,1)=i;
       end
       
     else
       for i=1:numel(CI)
           SI1(i,1)=i;
       end
       for i=1:numel(CI)
           for j=1:RL-1
               rs=(SI(ind1(i,j))-SI(ind1(i,j+1)))/SI(ind1(i,j));
               if rs <RS ...
                       | rs == RS
                  SI1(i,j+1)=ind1(i,j+1);
               else
                  SI1(i,j+1)=ind1(i,j);
               end
           end
       end
  
     end
   
       %creat rule
       y1=zeros(1,RL);
       p=ones(1,RL);
       for i=1:SI1(:,1)
           for j=1:RL
               y1(i)=y1(i)+a(SI1(i,j))*WF(SI1(i,j));
           end
           y1(i)=y1(i)+x(end-(5-i));
       end
       
       %aggregation
       for j=1:RL
        for i=1:SI1  
           p(j)=p(j)*SI(SI1(i,j));
        end
       end
       y=sum(p.*y1);
       

end