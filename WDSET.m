function [OutputScores,OutputLabels]=DSET(InputScores)
pos_prob=InputScores;
[no_sample,no_mass]=size(pos_prob);
neg_prob=1-pos_prob;
num_Sensor=no_mass;   
num_Object=4;  
num_Period=no_sample;   

Info=zeros(num_Sensor,num_Object,num_Period);

for i=1:num_Period
    for j=1:num_Sensor
     Info(j,:,i)=[pos_prob(i,j), neg_prob(i,j),0,0];
    end
end
 Info1=zeros(num_Period,num_Object); 
 k_threshold=0.5;
    for i=1:num_Period 
        Info1(i,:)=Info(1,:,i); 
        for j=1:num_Sensor-1 
            [Info1(i,:),k(i,:)]=DS_fusion(Info1(i,:),Info(j+1,:,i));
        end 
            if k(i)>k_threshold
                 innerProdm=zeros(no_mass,no_mass);
                     for s=1:no_mass
                         for t=1:no_mass
                             innerProdm(s,t)=Info(s,1,i)*Info(t,1,i)+Info(s,2,i)*Info(t,2,i);
                         end
                     end
                     dm=zeros(no_mass,no_mass);
                     for s=1:no_mass
                         for t=1:no_mass
                             dm(s,t)=(0.5*(innerProdm(s,s)+innerProdm(t,t)-2*innerProdm(s,t)))^0.5;
                         end
                     end
                     Simm=1-dm;
                     Supm=sum(Simm,2)-1;
                     Crdm=Supm./sum(Supm);
                    weighted_info=Info(:,:,i).*repmat(Crdm,1,num_Object);
                     weighted_info(:,3)=1-weighted_info(:,1)-weighted_info(:,2);
                     Info1(i,:)=weighted_info(1,:,1); 
                     for j=1:num_Sensor-1 
                        [Info1(i,:),kk(i,:)]=DS_fusion(Info1(i,:),weighted_info(j+1,:)); 
                     end
            end
    end 
for i=1:no_sample
    Result=Info1(i,:);
    [nx,mx]=size(Result); 
    if 1~=nx 
        disp('Result must be row vector'); 
        return; 
    end 
    ec1=0.05;
    ec2=0.1;
    [data,index]=sort(Result(1,1:mx-2));  
    if (Result(index(mx-2))-Result(index(mx-3))>ec1)  &&  (Result(mx-1)<ec2)  &&  (Result(index(mx-2))>Result(mx-1)) 
                if index(mx-2)==1
                    pred_label(i,1)=1;
                    DS_prob(i,1)=1-Result(2); 
                elseif index(mx-2)==2
                    pred_label(i,1)=-1;
                    DS_prob(i,1)=Result(1); 
                end
    else
                 DS_prob(i,1)=mean(Info(:,1,i));
                 if DS_prob(i,1)>=0.5
                    pred_label(i,1)=1;
                 else
                     pred_label(i,1)=-1;
                 end        
    end
end
OutputScores=DS_prob;
OutputLabels=pred_label;
end
