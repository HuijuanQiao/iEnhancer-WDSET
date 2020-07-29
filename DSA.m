clc
clear 
tt1=fastaread('D:\qiaohuijuan\test\enhancer\independent\independent dataset\enhancers.txt');
tt2=fastaread('D:\qiaohuijuan\test\enhancer\independent\independent dataset\nonhancers.txt');
tt=[tt1;tt2];
G=10;
n=length(tt);
feature2=zeros(n,G*36);

for i=1:n
    SSS{1,i}=tt(i).Sequence;
    Seq=SSS{1,i};
    
    %生成DNA样本的特征向量 
    value_na = cell(1,16);
    value_na{1} = 'AA'; value_na{2} = 'TC'; value_na{3} = 'AG'; value_na{4} = 'AT';
    value_na{5} = 'CA'; value_na{6} = 'CC'; value_na{7} = 'CG'; value_na{8} = 'CT';
    value_na{9} = 'GA'; value_na{10} = 'GC'; value_na{11} = 'GG'; value_na{12} = 'GT';
    value_na{13} = 'TA'; value_na{14} = 'TC'; value_na{15} = 'TG'; value_na{16} = 'TT';
    
    Pro = [7.65 8.93 7.08 9.07 6.38 8.04 6.23 7.08 8.56 9.53 8.04 8.93 6.23 8.56 6.38 7.65;
        2.26 3.03 2.03 3.83 1.78 1.65 2.00 2.03 1.93 2.61 1.65 3.03 1.20 1.93 1.78 2.26;
        1.69 1.32 1.46 1.03 1.07 1.43 1.08 1.46 1.32 1.20 1.43 1.32 0.72 1.32 1.07 1.69;
        0.026 0.036 0.031 0.033 0.016 0.026 0.014 0.031 0.025 0.025 0.026 0.036 0.017 0.025 0.016 0.026;
        0.020 0.023 0.019 0.022 0.017 0.019 0.016 0.019 0.020 0.026 0.019 0.023 0.016 0.020 0.017 0.020;
        0.038 0.038 0.037 0.036 0.025 0.042 0.026 0.037 0.038 0.036 0.042 0.038 0.018 0.038 0.025 0.038];
    
    for ii = 1:size(Pro,1)
        Pro(ii,:) = (Pro(ii,:)-min(Pro(ii,:)))/(max(Pro(ii,:))-min(Pro(ii,:)));
    end
    
    value_matrix = zeros(size(Pro,1),length(Seq)-1);
    for ii = 1:size(Pro,1)
        for jj = 1:16
            value_matrix(ii,strfind(Seq,value_na{jj})) = Pro(ii,jj);
        end
    end
    value_matrix = value_matrix';
   % G = 9;
    M = zeros(6*6,G);
    FF = zeros(6);
    for g = 1:G
        for s = 1:6
            for t = 1:6
            Csg = value_matrix(1:length(Seq)-1-g,s);
            Csg1= value_matrix(1+g:length(Seq)-1,s) ;
            Ctg= value_matrix(1:length(Seq)-1-g,t) ;
            Ctg1 = value_matrix((1+g):(length(Seq)-1),t);
            C1=value_matrix(:,s);
            C2=value_matrix(:,t);       
            t1=sum((Csg-Csg1).*(Ctg-Ctg1));
            l1=2*(length(Seq)-1-g);
            t2=sum(C1-mean(C1).*(C2-mean(C2)));
            l2=(length(Seq)-1-1);            
            FF(s,t) =(t1/l1)/(t2/l2);

            end
        end
        M(:,g) = FF(:);
    end
    feature2(i,:) = M(:)';
end
csvwrite('D:\qiaohuijuan\test\enhancer\independent\features\DSA360.csv',feature2);