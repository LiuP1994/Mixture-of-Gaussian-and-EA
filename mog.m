% Simple script to do EM for a mixture of Gaussians.
% -------------------------------------------------
%  based on code from  Rasmussen and Ghahramani
% (http://www.gatsby.ucl.ac.uk/~zoubin/course02/)

% Initialise parameters
%generate a data set x that contains only the F1 and F2
load('PB_data.mat')
for i = find(phno == 1)
    first_f1 = f1(i);
    first_f2 = f2(i);
end
phno1 = [first_f1, first_f2];

for i = find(phno == 2)
    second_f1 = f1(i);
    second_f2 = f2(i);
end
phno2 = [second_f1, second_f2];
x = [f1, f2, f1+f2];                    

[n D] = size(x);       
k = 6;                              

p = ones(1,k)/k;         

mu = x(ceil(n.*rand(1,k)),:)';     

s2 = zeros(D,D,k);                  

niter=100;                   

for i=1:k
  s2(:,:,i) = cov(x)./k;      
end

set(gcf,'Renderer','zbuffer');

clear Z;
try
  % run EM for niter iterations
  for t=1:niter,
    fprintf('t=%d\r',t);
    % Do the E-step:
    for i=1:k
      Z(:,i) = p(i)*det(s2(:,:,i))^(-0.5)*exp(-0.5*sum((x'-repmat(mu(:,i),1,n))'*inv(s2(:,:,i)).*(x'-repmat(mu(:,i),1,n))',2));
    end
    Z = Z./repmat(sum(Z,2),1,k);
    for i=1:k
      mu(:,i) = (x'*Z(:,i))./sum(Z(:,i));
      s2(:,:,i) = diag((x'-repmat(mu(:,i),1,n)).^2*Z(:,i)./sum(Z(:,i))); 
      p(i) = mean(Z(:,i));
    end
    
    clf
    hold on
    plot(x(:,1),x(:,2),'.');
    for i=1:k
      plot_gaussian(2*s2(:,:,i),mu(:,i),i,11);
    end
    
    drawnow;
  end
  
catch
  disp('Numerical Error in Loop - Possibly Singular Matrix');
end


%% Task3

for i = find(phno==1)
    F1 = f1(i);
    F2 = f2(i);
end
x_1=[F1,F2];

for i = find(phno==2)
    F1 = f1(i);
    F2 = f2(i);
end
x_2=[F1,F2];

y = zeros(152,2);
for i=1:length(y)
    y(i,1) = 1;
    y(i,2) = 2; 
end

output = zeros(152,2);

[n,D] = size(x_1);                    % number of observations (n) and dimension (D)
k = 3;                              % number of components

clear Z;
%probabilities for model 1 phenome 1
zsum_1 = zeros(152,2);
zsum_2 = zeros(152,2);
  for i=1:k
      Z(:,i) = p_1(i)*det(s2_1(:,:,i))^(-0.5)*exp(-0.5*sum((x_1'-repmat(mu_1(:,i),1,n))'*inv(s2_1(:,:,i)).*(x_1'-repmat(mu_1(:,i),1,n))',2));
      zsum_1(:,1) = zsum_1(:,1) + Z(:,i);
  end
  for i=1:k
      Z(:,i) = p_1(i)*det(s2_1(:,:,i))^(-0.5)*exp(-0.5*sum((x_2'-repmat(mu_1(:,i),1,n))'*inv(s2_1(:,:,i)).*(x_2'-repmat(mu_1(:,i),1,n))',2));
      zsum_2(:,1) = zsum_2(:,1) + Z(:,i);
  end
  for i=1:k
      Z(:,i) = p_2(i)*det(s2_2(:,:,i))^(-0.5)*exp(-0.5*sum((x_1'-repmat(mu_2(:,i),1,n))'*inv(s2_2(:,:,i)).*(x_1'-repmat(mu_2(:,i),1,n))',2));
      zsum_1(:,2) = zsum_1(:,2) + Z(:,i);
  end
  for i=1:k
      Z(:,i) = p_2(i)*det(s2_2(:,:,i))^(-0.5)*exp(-0.5*sum((x_2'-repmat(mu_2(:,i),1,n))'*inv(s2_2(:,:,i)).*(x_2'-repmat(mu_2(:,i),1,n))',2));
      zsum_2(:,2) = zsum_2(:,2) + Z(:,i);
  end

for i=1:length(x)
    if zsum_1(i,1) > zsum_1(i,2)
        output(i,1) = 1;
    else
        output(i,1)= 2;
    end
end
for i=1:length(x)
    if zsum_2(i,1) > zsum_2(i,2)
        output(i,2) = 1;
    else
        output(i,2)= 2;
    end
end

miss_classification1 = 1-sum(output(:,1)==1)/152;
miss_classification2 = 1-sum(output(:,2)==2)/152;
miss_classification = 1-(sum(output(:,1)==1)+sum(output(:,2)==2))/304;
disp("miss_classification1="+miss_classification1)
disp("miss_classification2="+miss_classification2)
disp("miss_classification="+miss_classification)


%% Task 4

for i = find(phno==1)
    F1 = f1(i);
    F2 = f2(i);
end

for i = find(phno==2)
    F3 = f1(i);
    F4 = f2(i);
end
k = 6; 
x1=[F1,F3];
x2=[F2,F4];

for i = min(x1):max(x1)
    for j = min(x2):max(x2)
        sum1 = 0;
        sum2 = 0;
        x = [i,j];
        for t=1:k
            n1(:,t) = p_1(t)*det(s2_1(:,:,t))^(-0.5)*exp(-0.5*sum((x'-repmat(mu_1(:,t),1))'*inv(s2_1(:,:,t)).*(x'-repmat(mu_1(:,t),1))',2)); 
            sum1 = sum1 + n1(:,t);
 
            n2(:,t) = p_2(t)*det(s2_2(:,:,t))^(-0.5)*exp(-0.5*sum((x'-repmat(mu_2(:,t),1))'*inv(s2_2(:,:,t)).*(x'-repmat(mu_2(:,t),1))',2)); 
            sum2 = sum2 + n2(:,t);
        end
        if sum1 > sum2
            values = 1;
        else
            values = 2;
        end
        M(i,j) = values;
    end
end
imagesc(M')

