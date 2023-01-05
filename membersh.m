function y= membersh(x,A,mo,spread,Data)

    global CI
    %w
    k=1;
    w=zeros(numel(CI),CI(1));
    for j=1:numel(CI)
        for i=1:CI(1)
            w(j,i)=x(k);
            k=k+1;
           
        end
    end
    
    %a
    a=zeros(numel(CI),CI(1));
    for j=1:numel(CI)
        for i=1:CI(1)
            a(j,i)=x(k);
            k=k+1;
        end
    end
    
    %RL,RS
    RL=round(x(end-1));
    RS=x(end);
    
    % MDOV matrix
    MODV=zeros(numel(CI),CI(1));
    for d=1:size(Data,1)
      for i=1:numel(CI)
        for j=1:CI(1)
            MODV(i,j)=w(i,j)*gaussmf(Data(d,i),[spread(j) mo(j)]);
        end
      end
      
    %SI
    for i=1:numel(CI)
        SI(i,1)=MODV(i,1)+MODV(i,2);
        SI(i,CI)=MODV(i,CI)+MODV(i,CI-1);
    end
    for j=1:numel(CI)
       for i=2:CI-1
        SI(j,i)=MODV(j,i-1)+MODV(j,i)+MODV(j,i+1);
       end
    end
    
    
    
    [MODV ind]=sort(MODV,2,'descend');% sort from big to short
      
    % calculate sorted index
    ind1=zeros(numel(CI),CI(1));
    for d1=1:numel(CI)
          for i1=1:CI(1)
              for k=1:CI(1)
                  if ind(d1,i1)==k
                      ind1(d1,i1)=(k-1)*numel(CI)+d1;
                  end
              end
          end
    end
     
    % creat rules
    y(d)=rules(x,ind1,MODV,a,RS,3,SI);
      
      
    end
    
 
 end