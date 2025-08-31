%   A reliable and adaptive spatiotemporal data fusion (RASDF) method

%   Using one pairs of fine and coarse images;
%   Inputing images of envi format, resampling the coarse images to make it consistent with the size of fine image.


%   The author would like to thank Prof. Xiaolin Zhu, Dr. Xiaodong Li, Dr. Yue Sun, et al 
%   for providing the source code of FSDAF, SFSDAF and OBSTFM, which were referenced in this program 

%   RASDF version 1.0  02/10/2021
%   RASDF version 1.1  18/01/2022 Possible bugs in change detection and edge detection modification 

clear all
clc

%%%%%%%%%%%%%%%%   Parameters  %%%%%%%%%%%%%%%%%
%Parameters needed setting
DN_min = 0;  % Minimum DN value
DN_max = 10000;   % Maximum DN value; if using reflectance, use 1 as DN_max
Factor=20;  % Resolution ratio of coarse image to fine image

%The following parameters are recommended as default values:
options=[2,80,1e-6,true];   % Parameters of FCM in Global unmixing
options2=[2,50,1e-8,true];  % Parameters of FCM in Local unmixing
min_class=4;    % Minimum number of classification in Global unmixing step
max_class=7;    % Maximum number of classification in Global unmixing step
min_class2=2;   % Minimum number of classification in Local unmixing step
max_class2=4;   % Maximum number of classification in Local unmixing step
w_min=20; % Minimum window size in filtering and residual compensation steps
w_max=80; % Maximum window size in filtering and residual compensation steps
ww2=2;  % window size in Local unmixing step
min_num=max(9,max_class2);

%%%%%%%%%%%%%%%%   The path of the images  %%%%%%%%%%%%%%%%%
tif='.tif';
tiff='.tiff';
filename1 = 'E:\已发表论文或不用的实验材料\dongtinlakecode\20240707flood\forfusion\LC08_124040_20240517_cut.tif';
[~,~,format]=fileparts(filename1);
if strcmp(format,tif) || strcmp(format,tiff)
    FRT1 = imread(filename1);
else
    FRT1 = enviread(filename1);
end
FRT1(FRT1<DN_min)=DN_min;

filename2 = 'E:\dongtinlakecode\20240707flood\forfusion\MCD43A42024_05_17_cut.tif';
[~,~,format]=fileparts(filename2);
if strcmp(format,tif) || strcmp(format,tiff)
    CRT1 = imread(filename2);
else
    CRT1 = enviread(filename2);
end
CRT1(CRT1<DN_min)=DN_min;

filename3 = 'E:\dongtinlakecode\20240707flood\forfusion\MOD02HK2024_07_07_cut.tif';
[~,~,format]=fileparts(filename3);
if strcmp(format,tif) || strcmp(format,tiff)
    CRT2 = imread(filename3);
else
    CRT2 = enviread(filename3);
end
CRT2(CRT2<DN_min)=DN_min;

outputname = 'E:\dongtinlakecode\20240707flood\forfusion\RASDF2024_07_07';


[xH,yH,bands]=size(FRT1);

CRT1=double(CRT1);
CRT2=double(CRT2);
FRT1=double(FRT1);

%%%%%%%%%%%%%%%%   Resampling  %%%%%%%%%%%%%%%%%
xL=xH/Factor;
yL=yH/Factor;

CRT1C=zeros(xL,yL,bands);
CRT2C=zeros(xL,yL,bands);
FRT1C=zeros(xL,yL,bands);
for i=1:xL
    for j=1:yL
        tmp1 = CRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        CRT1C(i,j,:)=sum(sum(tmp1))/Factor^2;
        tmp2 = CRT2( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        CRT2C(i,j,:)=sum(sum(tmp2))/Factor^2;
        tmp3 = FRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        FRT1C(i,j,:)=sum(sum(tmp3))/Factor^2;
    end
end

nearest_FRT1C = (imresize(FRT1C,Factor,'nearest'));
cubic_CR = (imresize(CRT2C-CRT1C,Factor,'cubic'));
cubic_CRT1=(imresize(CRT1C,Factor,'cubic'));

clear tmp1 tmp2 tmp3


%%%%%%%%%%%%%%%%   Obtain RI  %%%%%%%%%%%%%%%%%
RI=zeros(xH,yH,bands);

X=reshape(CRT1C,xL*yL,bands);
Y=reshape(FRT1C,xL*yL,bands);

for bb=1:bands
    x=zeros(xL*yL,1);
    y=zeros(xL*yL,1);
    x(:)=X(:,bb);
    y(:)=Y(:,bb);
    xx=[ones(length(y),1),x];
    [BB,Bint,R]=regress(y,xx);
    a(bb)=BB(2);
    b(bb)=BB(1);
    CRT1CN(:,:,bb)=a(bb).*CRT1(:,:,bb)+b(bb);
    CRT2CN(:,:,bb)=a(bb).*CRT2(:,:,bb)+b(bb);
end


Fd=CRT1CN-nearest_FRT1C;
Td=CRT2-CRT1;

for b=1:bands
    mean_Fd=mean2(abs(Fd(:,:,b)));
    mean_Td=mean2(abs(Td(:,:,b)));
    std_Fd=std2(Fd(:,:,b));
    
    Q_m(b)=2*mean_Td./mean_Fd;
    Q_m(Q_m<2)=2;
    for i=1:xH
        for j=1:yH
            if abs(Fd(i,j,b)) < Q_m(b)*std_Fd
                RI(i,j,b)=1-abs(Fd(i,j,b))/(Q_m(b)*std_Fd);
            end
        end
    end
    
end
RI(RI<0.1)=0.1;
RIC=zeros(xL,yL,bands);
for i=1:xL
    for j=1:yL
        tmp = RI( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: ); 
        RIC(i,j,:)=sum(sum(tmp))/Factor^2;
    end
end


w_s=round(w_max-mean(mean(mean(RIC)))*(w_max-w_min));
similar_pixel=w_s;

%%%%%%%%%%%%%%%  Change detection for Global unmix  %%%%%%%%%%%%%%
%   Gaussian
disp('Change detection for Global unmix')
Change_detection_unmix=zeros(xL,yL);
Change_detection_unmix2=zeros(xL,yL,bands)+1;
Q_threshold=zeros(bands,1);
Q_mean=zeros(bands,1);

D_C=CRT2C-CRT1C;

for b=1:bands
    D_C_b=D_C(:,:,b);
    D_C_b_s=reshape(D_C_b,1,xL*yL);
    
    Q_threshold(b)=2*std2(D_C_b(:,:));
    Q_mean(b)=mean2(D_C_b(:,:));

    min_allow(b)=Q_mean(b)-Q_threshold(b);
    max_allow(b)=Q_mean(b)+Q_threshold(b);
    for i=1:xL
        for j=1:yL
            if ((D_C(i,j)>max_allow(b))||(D_C(i,j)<min_allow(b)))
                Change_detection_unmix(i,j)=1;
                Change_detection_unmix2(i,j,:)=0;
            end
        end
    end
end

%OTSU
% for b=1:bands
%     D_C_b=D_C(:,:,b);
%     D_C_b_s=reshape(D_C_b,1,xL*yL);
%     
%      if (mean2(Change_detection_unmix)>=0.75)&&(b==bands-1) % if input 
%         Change_detection_unmix=zeros(xL,yL);
%         Change_detection_unmix2=zeros(xL,yL,bands)+1; %因为后面要乘它，表示变化区域为0，只在非变化区域判断
%         D_C_pos=D_C_b_s(find(D_C_b_s>=0));
%         D_C_neg=D_C_b_s(find(D_C_b_s<0));
%         
%         I1=abs(D_C_pos)/DN_max;
%         I2=abs(D_C_neg)/DN_max;
%         
%         max_allow(b)=DN_max*graythresh(I1);
%         min_allow(b)=-DN_max*graythresh(I2);
%         
%         for i=1:xL
%             for j=1:yL
%                 if ((D_C_b(i,j)>max_allow(b))||(D_C_b(i,j)<min_allow(b)))
%                     Change_detection_unmix(i,j)=1;
%                     Change_detection_unmix2(i,j,:)=0;
%                 end
%             end
%         end
%         
%      end
% end
Change_detection_unmix_f=zeros(xH,yH);

for b=1:bands
    for i=1:xL
        for j=1:yL    
            Change_detection_unmix_f( (i-1)*Factor+1:i*Factor,(j-1)*Factor+1:j*Factor)=Change_detection_unmix(i,j);
        end
    end
end

%%%%%%%%%%%%%%%%%   Clustering  %%%%%%%%%%%%%%%%
disp('Cluster T1 FR image')
dif1=min_class-1;
FR_T1_V = zeros(xH*yH,bands);
F2TP_class=zeros(xH,yH,bands,max_class-dif1);
CR_fraction_T1_class=zeros(xL,yL,max_class,max_class-dif1);
FR_fraction_T1_class=zeros(xH,yH,max_class,max_class-dif1);

for i=1:bands
    FR_T1_V(:,i) = reshape(FRT1(:,:,i),xH*yH,1);
end  
for class=min_class:max_class
    CR_fraction_T1=zeros(xL,yL,class);
    [~, result] = fcm(FR_T1_V,class,options);
    map_T1=reshape(result',xH,yH,class);
    [~,result3]=max(result,[],1);
    clear result 
    FR_fraction_T1=map_T1;
     for c=1:class
         for i=1:xL
             for j=1:yL
                 tmp = sum(sum(  FR_fraction_T1( (i-1)*Factor+1:i*Factor,(j-1)*Factor+1:j*Factor,c) ))/Factor^2;
                 CR_fraction_T1(i,j,c)=tmp;
                 CR_fraction_T1_class(i,j,c,class-dif1)=tmp;
                 FR_fraction_T1_class(i,j,c,class-dif1)=FR_fraction_T1(i,j,c);
             end
         end
     end
    clear tmp 
    
    %%%%%%%%%%%%%%%%%%%%%%%  Global unmix  %%%%%%%%%%%%%%%%%%%%%%%
    disp('Estimate T1 to T2 FR endmember change')
    FR_endmember_change=zeros(bands,class);
    for b=1:bands
        CRT2Cb=CRT2C(:,:,b);
        CRT1Cb=CRT1C(:,:,b);
        min_allowb=min_allow(b);
        max_allowb=max_allow(b);
        RICb=RIC(:,:,b);
        RICb=reshape(RICb,1,xL*yL);
        [RI_sort,~]=sort(RICb,'descend');
        Q_RI=RI_sort(round(0.3*xL*yL));
        
        
        mmm=1;
        nnn=1;
        mm=1;
        nn=1;
        p=1;
        pp=1;
        for i=1:xL
            for j=1:yL                
                if (RIC(i,j,b)<=Q_RI)&&(RIC(i,j,b)<0.95)
                    mmm(pp)=i;
                    nnn(pp)=j;
                    pp=pp+1;
                end
                if Change_detection_unmix(i,j)==1
                    mm(p)=i;
                    nn(p)=j;
                    p=p+1;
                end
            end
        end
        
        tmp=Estimate_endmember_change_each_band(CR_fraction_T1 ,CRT2Cb ,CRT1Cb ,xL,yL,class,1,mm,nn,mmm,nnn,min_allowb,max_allowb,p,pp);
        FR_endmember_change(b,:)=tmp(:);
    end

    FchangeC=zeros(xL,yL,bands);
    temp=zeros(xH,yH,bands);
    for b=1:bands
        for c=1:class
            temp(:,:,b) = temp(:,:,b) + FR_endmember_change(b,c) * FR_fraction_T1(:,:,c) ;
        end
    end
    F2TP=FRT1+temp; 
    for i=1:xL
        for j=1:yL
            tmp = F2TP( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: )-FRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:); 
            FchangeC(i,j,:)=sum(sum(tmp))/Factor^2;
        end
    end
    
    F2TP_class(:,:,:,class-dif1)=F2TP(:,:,:);    
    residuals_c_sum=reshape(sum(sum(abs((CRT2C-CRT1C)-FchangeC)./abs((CRT2C-CRT1C)+0.001).*RIC.*Change_detection_unmix2)),1,bands);
    residuals_c_sum_class(class-dif1)=sum(residuals_c_sum);

    clear temp
end

   [~,index1]=min(residuals_c_sum_class(:));
    for i=1:xL
        for j=1:yL
            F2TP((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)=F2TP_class((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:,index1);
        end
    end

CR_fraction_T1_best=reshape(CR_fraction_T1_class(:,:,1:index1+min_class-1,index1),xL,yL,index1+min_class-1);
FR_fraction_T1_best=reshape(FR_fraction_T1_class(:,:,1:index1+min_class-1,index1),xH,yH,index1+min_class-1);

XX = sprintf('The number of classifiation in Global unmixing:%d',index1+min_class-1);
disp(XX)

clear  F2TP_class index FchangeC FR_fraction_T1

%%%%%%%%%%%%%%%%  Local unmixing  %%%%%%%%%%%%%%%%%%%%%

disp('Local unmixing')
dif2=min_class2-1;  
FRT2s=FRT1+cubic_CR; 
FRT2s2=nearest_FRT1C+cubic_CR; 
FchangeC=zeros(xL,yL,bands);

for i=1:xL
    for j=1:yL
        tmp = F2TP( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: )-FRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:); 
        FchangeC(i,j,:)=sum(sum(tmp))/Factor^2;
    end
end
CCR=((CRT2C-CRT1C)-FchangeC);

change_detection_c=zeros(xL,yL);
for b=1:bands 
    CCR_std(b)=std2(CCR(:,:,b));
    CCR_mean(b)=mean2(CCR(:,:,b));
    for i=1:xL
        for j=1:yL
            if (CCR(i,j,b)>CCR_mean(b)+CCR_std(b)*2)||(CCR(i,j,b)<CCR_mean(b)-CCR_std(b)*2)
                change_detection_c(i,j)=1; 
            end
        end
    end
end


RIC_mb=mean(RIC,3);

C_change=CRT2C-CRT1C;
F2TP2=F2TP;
F2TP3=F2TP;
number2=zeros(xL,yL);
class_index=zeros(xL,yL);

for i=1:xL
    for j=1:yL
         if (change_detection_c(i,j)==1)&&((mean2(RIC(i,j,:))>0.8))
            ai = max(1,i-ww2);
            bi = min(xL,i+ww2);
            aj = max(1,j-ww2);
            bj = min(yL,j+ww2);
            SIC_in = reshape(RIC(ai:bi,aj:bj,:),(bi-ai+1)*(bj-aj+1),bands);
            SIC_mbin=reshape(RIC_mb(ai:bi,aj:bj),(bi-ai+1)*(bj-aj+1),1);
            CRT2C_in = reshape(CRT2C(ai:bi,aj:bj,:),(bi-ai+1)*(bj-aj+1),bands);
            CRT1C_in = reshape(CRT1C(ai:bi,aj:bj,:),(bi-ai+1)*(bj-aj+1),bands);
            change_detection_c_in=change_detection_c(ai:bi,aj:bj);
            change_detection_c_in=change_detection_c_in(:);
            spediff1=sum(abs(C_change(ai:bi,aj:bj,:)-C_change(i,j,:))./C_change(i,j,:),3);
            spediff1=spediff1(:);
            Q_C=std(spediff1)/2;
            starfm_similar_index1=find(spediff1  < Q_C & change_detection_c_in==1 & SIC_mbin>=0.8);
            num_starfm_similar1=length(starfm_similar_index1(:));
            row= mod(starfm_similar_index1,(bi-ai+1))+ai-1; 
            col= fix(starfm_similar_index1/(bi-ai+1))+aj; 
            number2(i,j)=num_starfm_similar1;
            for row_i=1:num_starfm_similar1  
                if row(row_i)==ai-1
                    row(row_i)=bi;
                    col(row_i)=col(row_i)-1;
                end            
            end
            if num_starfm_similar1>=min_num    % Case 1 in local unmixing
                CRT2C_in2 = CRT2C_in(starfm_similar_index1,:);
                CRT1C_in2 = CRT1C_in(starfm_similar_index1,:);
                SIC_in2=SIC_in(starfm_similar_index1,:);
               
                FR_T1_V1 = zeros(num_starfm_similar1*Factor,Factor,bands);
                FR_T1_V2=zeros(num_starfm_similar1*Factor*Factor,bands);
                for o=1:num_starfm_similar1
                    FR_T1_V1((o-1)*Factor+1:o*Factor,1:Factor,:)=FRT1((row(o)-1)*Factor+1:row(o)*Factor, (col(o)-1)*Factor+1:col(o)*Factor,:);
                    if row(o)==i && col(o)==j
                        FRT1_target=FRT1((row(o)-1)*Factor+1:row(o)*Factor, (col(o)-1)*Factor+1:col(o)*Factor,:);
                        o_target=o;
                    end
                end                
                new_targetw=zeros(num_starfm_similar1,bands)+(0.5)/(num_starfm_similar1);
                
                 for bb=1:bands
                      FR_T1_V2(:,bb) = reshape(FR_T1_V1(:,:,bb),num_starfm_similar1*Factor*Factor,1);
                 end  
                 
                for class=min_class2:max_class2
                    [~, result2] = fcm(FR_T1_V2,class,options2);
                    map_T1_in=reshape(result2',num_starfm_similar1*Factor,Factor,class);                   
                    clear result2
                    CR_fraction_T1=zeros(num_starfm_similar1,class);
                    FR_fraction_T1=map_T1_in;
                     for c=1:class
                         for ii=1:num_starfm_similar1                             
                             tmp = sum(sum(  FR_fraction_T1( (ii-1)*Factor+1:ii*Factor,1:Factor,c) ))/Factor^2;
                             CR_fraction_T1(ii,1,c)=tmp;                             
                         end
                     end
                     FR_endmember_change = unmix2(CR_fraction_T1 ,CRT2C_in2 ,CRT1C_in2 ,num_starfm_similar1,1,class,bands);
                     temp=zeros(num_starfm_similar1*Factor,Factor,bands);
                     for b=1:bands
                         for c=1:class
                             temp(:,:,b) = temp(:,:,b) + FR_endmember_change(b,c) * FR_fraction_T1(:,:,c) ;
                         end
                     end
                     F2TP_in=FRT1_target+temp((o_target-1)*Factor+1:o_target*Factor, 1:Factor,:); 
                      temp_in=zeros(num_starfm_similar1,bands);
                      for iii=1:num_starfm_similar1
                          tmp = temp( (iii-1)*Factor+1:iii*Factor, 1:Factor,: ); 
                          temp_in(iii,:)=sum(sum(tmp))/Factor^2;
                      end
                      CCR_in=sum(sum(abs((CRT2C_in2-CRT1C_in2)-temp_in).*SIC_in2.*new_targetw));              
                      C_centrol(class-dif2)=CCR_in;
                      F2TP_class(1:Factor,1:Factor,:,class-dif2)=F2TP_in;
                end

                CR_fraction_T1_in=reshape(CR_fraction_T1_best(ai:bi,aj:bj,:),(bi-ai+1)*(bj-aj+1),index1+min_class-1);            
                CR_fraction_T1=CR_fraction_T1_in(starfm_similar_index1,:);
                CR_fraction_T1=reshape(CR_fraction_T1,num_starfm_similar1,1,index1+min_class-1);
                for c=1:index1+min_class-1
                    for ii=1:num_starfm_similar1
                        FR_fraction_T1((ii-1)*Factor+1:ii*Factor,1:Factor,c)=FR_fraction_T1_best((row(o)-1)*Factor+1:row(o)*Factor, (col(o)-1)*Factor+1:col(o)*Factor,c);
                    end
                end
                FR_endmember_change = unmix2(CR_fraction_T1 ,CRT2C_in2 ,CRT1C_in2 ,num_starfm_similar1,1,index1+min_class-1,bands);
                temp=zeros(num_starfm_similar1*Factor,Factor,bands);
                for b=1:bands
                    for c=1:index1+min_class-1
                        temp(:,:,b) = temp(:,:,b) + FR_endmember_change(b,c) * FR_fraction_T1(:,:,c) ;
                    end
                end
                F2TP_in=FRT1_target+temp((o_target-1)*Factor+1:o_target*Factor, 1:Factor,:); 
                temp_in=zeros(num_starfm_similar1,bands);
                for iii=1:num_starfm_similar1
                    tmp = temp( (iii-1)*Factor+1:iii*Factor, 1:Factor,: ); 
                    temp_in(iii,:)=sum(sum(tmp))/Factor^2;
                end
                CCR_in=sum(sum(abs((CRT2C_in2-CRT1C_in2)-temp_in).*SIC_in2.*new_targetw));                 
                C_centrol(class-dif2+1)=CCR_in;
                F2TP_class(1:Factor,1:Factor,:,class-dif2+1)=F2TP_in;
                [~,index]=min(abs(C_centrol(1:class-dif2+1)));
                F2TP2((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)=F2TP_class(1:Factor,1:Factor,:,index);
                class_index(i,j)=index+min_class2-1;

                ccc=sum(abs((CRT2C(i,j,:)-CRT1C(i,j,:)) - (sum(sum(F2TP2((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:))))/(Factor^2)+FRT1C(i,j,:)));
                ddd=sum(abs(CRT2C(i,j,:)-CRT1C(i,j,:)));
                if (ccc<ddd)&&((((ddd-ccc)/ddd)>=0.7)||(((ddd-ccc)/ddd)>=mean(RIC(i,j,:)))) %Here, the threshold of 0.7 is to deal with the case when feeding simulated images (i.e., using downsampling images as coarse images), in this case, mean (RI) is close to 1.
                    F2TP3((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)=((ddd-ccc)/ddd)*(F2TP2((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)-FRT1((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:))+FRT1((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:);                    
                else
                    F2TP3((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)=F2TP((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:);
                end
        elseif num_starfm_similar1<min_num   % Case 2 in local unmixing
            ai = max(1,i-1);
            bi = min(xL,i+1);
            aj = max(1,j-1);
            bj = min(yL,j+1);
            CRT2C_in = CRT2C(ai:bi,aj:bj,:);
            CRT1C_in = CRT1C(ai:bi,aj:bj,:);
            RIC_in3 = RIC(ai:bi,aj:bj,:);
            new_targetw=zeros(bi-ai+1,bj-aj+1,bands)+(0.5)/((bi-ai+1)*(bj-aj+1));
            FRT1_in=FRT1((ai-1)*Factor+1:bi*Factor, (aj-1)*Factor+1:bj*Factor,:);
            FR_T1_V = zeros(Factor*(bi-ai+1)*Factor*(bj-aj+1),bands);
            for b=1:bands
                 FR_T1_V(:,b) = reshape(FRT1_in(:,:,b),Factor*(bi-ai+1)*Factor*(bj-aj+1),1);
            end  
            for class=min_class2:max_class2

                [~, result2] = fcm(FR_T1_V,class,options2);
                map_T1_in=reshape(result2',Factor*(bi-ai+1),Factor*(bj-aj+1),class);
                clear result2
                CR_fraction_T1=zeros(bi-ai+1,bj-aj+1,class);
                FR_fraction_T1=map_T1_in;

                 for c=1:class
                     for ii=1:bi-ai+1
                         for jj=1:bj-aj+1
                             tmp = sum(sum(  FR_fraction_T1( (ii-1)*Factor+1:ii*Factor,(jj-1)*Factor+1:jj*Factor,c) ))/Factor^2;
                             CR_fraction_T1(ii,jj,c)=tmp;
                         end
                     end
                 end
                 FR_endmember_change = unmix(CR_fraction_T1 ,CRT2C_in ,CRT1C_in ,bi-ai+1,bj-aj+1,class,bands);
                 temp=zeros(Factor*(bi-ai+1),Factor*(bj-aj+1),bands);
                 for b=1:bands
                     for c=1:class
                         temp(:,:,b) = temp(:,:,b) + FR_endmember_change(b,c) * FR_fraction_T1(:,:,c) ;
                     end
                 end
                 F2TP_in=FRT1_in((i-ai+1-1)*Factor+1:(i-ai+1)*Factor, (j-aj+1-1)*Factor+1:(j-aj+1)*Factor,:)+temp((i-ai+1-1)*Factor+1:(i-ai+1)*Factor, (j-aj+1-1)*Factor+1:(j-aj+1)*Factor,:); 
                
                  temp_in=zeros((bi-ai+1),(bj-aj+1),bands);
                  for iii=1:(bi-ai+1)
                      for jjj=1:(bj-aj+1)
                          tmp = temp( (iii-1)*Factor+1:iii*Factor, (jjj-1)*Factor+1:jjj*Factor,: ); 
                          temp_in(iii,jjj,:)=sum(sum(tmp))/Factor^2;
                      end
                  end
                  CCR_in=sum(sum(sum(abs((CRT2C_in-CRT1C_in)-temp_in).*RIC_in3.*new_targetw)));                 
                  C_centrol(class-dif2)=CCR_in;
                  F2TP_class(1:Factor,1:Factor,:,class-dif2)=F2TP_in;
            end

            CR_fraction_T1=CR_fraction_T1_best(ai:bi,aj:bj,:);
            FR_fraction_T1=FR_fraction_T1_best((ai-1)*Factor+1:bi*Factor, (aj-1)*Factor+1:bj*Factor,:);
            FR_endmember_change = unmix(CR_fraction_T1 ,CRT2C_in ,CRT1C_in ,bi-ai+1,bj-aj+1,index1+min_class-1,bands);
                 temp=zeros(Factor*(bi-ai+1),Factor*(bj-aj+1),bands);
                 for b=1:bands
                     for c=1:index1+min_class-1
                         temp(:,:,b) = temp(:,:,b) + FR_endmember_change(b,c) * FR_fraction_T1(:,:,c) ;
                     end
                 end
                 F2TP_in=FRT1_in((i-ai+1-1)*Factor+1:(i-ai+1)*Factor, (j-aj+1-1)*Factor+1:(j-aj+1)*Factor,:)+temp((i-ai+1-1)*Factor+1:(i-ai+1)*Factor, (j-aj+1-1)*Factor+1:(j-aj+1)*Factor,:); 

                  temp_in=zeros((bi-ai+1),(bj-aj+1),bands);
                  for iii=1:(bi-ai+1)
                      for jjj=1:(bj-aj+1)
                          tmp = temp( (iii-1)*Factor+1:iii*Factor, (jjj-1)*Factor+1:jjj*Factor,: ); 
                          temp_in(iii,jjj,:)=sum(sum(tmp))/Factor^2;
                      end
                  end

                  CCR_in=sum(sum(sum(abs((CRT2C_in-CRT1C_in)-temp_in).*RIC_in3.*new_targetw)));
                 
                 C_centrol(class-dif2+1)=CCR_in;
                 F2TP_class(1:Factor,1:Factor,:,class-dif2+1)=F2TP_in;

                 [~,index]=min(abs(C_centrol(1:class-dif2+1)));
                 F2TP2((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)=F2TP_class(1:Factor,1:Factor,:,index);

                 class_index(i,j)=index+min_class2-1;
            
                 ccc=sum(abs((CRT2C(i,j,:)-CRT1C(i,j,:)) - (sum(sum(F2TP2((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:))))/(Factor^2)+FRT1C(i,j,:)));
                 ddd=sum(abs(CRT2C(i,j,:)-CRT1C(i,j,:)));   

                if (ccc<ddd)&&((((ddd-ccc)/ddd)>=0.7)||(((ddd-ccc)/ddd)>=mean(RIC(i,j,:)))) %Here, the threshold of 0.7 is to deal with the case when feeding simulated images (i.e., using downsampling images as coarse images), in this case, mean (RI) is close to 1
                    F2TP3((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)=((ddd-ccc)/ddd)*(F2TP2((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)-FRT1((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:))+FRT1((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:);
                else
                    F2TP3((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:)=F2TP((i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:);
                end

            end
        end
    end
end
clear F2TP_in F2TP_class temp CRT12C_in CRT2C_in F12TP_in temp_in FR_endmember_change FR_fraction_T1

F2TP2=F2TP3;


%%%%%%%%%%%%%%%%% Filter  %%%%%%%%%%%%%%%%%

disp('Filter')

B=zeros(xH,yH,similar_pixel);
Similar_index=zeros(xH,yH,similar_pixel);
endmember_change=F2TP2-FRT1;
prediction=zeros(xH,yH,bands);
mean_FRT1_allb=mean(mean(FRT1,1),2);
mean_FRT1_allb(:)=mean_FRT1_allb;
FRT1b=zeros(xH,yH,bands);
number=zeros(xH,yH);
for b=1:bands
    FRT1b(:,:,b)=FRT1(:,:,b)/mean_FRT1_allb(b);
end
FRT1_allb=sum(FRT1b,3);
Q_similar=std2(FRT1_allb(:,:))*2/(index1+min_class);

for i=1:xH
    for j=1:yH
        w=w_s;
        ai=max(1,i-w);
        bi=min(xH,i+w);
        aj=max(1,j-w);
        bj=min(yH,j+w);
        spediff=sum(abs(FRT1(ai:bi,aj:bj,:)-FRT1(i,j,:))./(FRT1(i,j,:)+0.01),3);                    
        spediff=spediff(:);
        starfm_similar_index=find(spediff  < Q_similar ); %index in window
        num_starfm_similar=length(starfm_similar_index(:));
        num=min(num_starfm_similar,similar_pixel);
        number(i,j)=num;
        [~,index]=sort(spediff);  
        index2=sort(index(1:similar_pixel)); %index in window
        row= mod(index2,(bi-ai+1))+ai-1; 
        col= fix(index2/(bi-ai+1))+aj; 
        for row_i=1:similar_pixel  
            if row(row_i)==ai-1
               row(row_i)=bi;
               col(row_i)=col(row_i)-1;
             end            
         end

         locadiff=zeros(1,similar_pixel)+10000000;
         for z=1:num
             locadiff(z)=(((row(z)-i)^2+(col(z)-j)^2)^0.5);
         end

          cijk=1.0./(1+locadiff/w);
          B(i,j,1:similar_pixel)=cijk/sum(cijk);
          B2(:)=B(i,j,:);
          
          for b=1:bands
              endmember_change_in=endmember_change(ai:bi,aj:bj,b);
              endmember_change_inx(:)=endmember_change_in(index2(:));
              prediction(i,j,b)=B2*endmember_change_inx'+FRT1(i,j,b);
          end  
           %Change 2D coordinates to 1D coordinates for easy storage and recall
           one_dis_index=(col-1)*xH+row;           
           %Save similar pixel locations
           for one=1:similar_pixel
                Similar_index(i,j,one)=one_dis_index(one);                
           end             
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%  Residual distribution  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RI = (imresize(RIC,Factor,'cubic')); %nearest

disp('Residual distribution')
FchangeC=zeros(xL,yL,bands);

for i=1:xL
    for j=1:yL
        tmp = prediction( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,: )-FRT1( (i-1)*Factor+1:i*Factor, (j-1)*Factor+1:j*Factor,:); 
        FchangeC(i,j,:)=sum(sum(tmp))/Factor^2;
    end
end
residuals_first=(CRT2C-CRT1C)-FchangeC;
residuals_f = imresize(residuals_first,Factor,'cubic');
mean_residuals_first=(sum(sum(abs(residuals_first))))/(xL*yL);
prediction2=zeros(xH,yH,bands);
for b=1:bands
    residuals_f_b=residuals_f(:,:,b);
    SI_b=RI(:,:,b);
    for i=1:xH
        for j=1:yH
            Gain=0;
            B3=reshape(B(i,j,:),1,similar_pixel);
            index=reshape(Similar_index(i,j,:),1,similar_pixel);
            sum_SI_b=sum(SI_b(index(1:number(i,j))));
            for one=1:number(i,j)
                SII=SI_b(index(one))/sum_SI_b;
                Gain=Gain+(B3(one)*residuals_f_b(index(one))*SII*number(i,j));
            end
            prediction2(i,j,b)=Gain+prediction(i,j,b);
        end
    end
end
first_residual=prediction2-prediction;
    
%Eliminate extreme values
for b=1:bands
    for i=1:xH
        for j=1:yH
            if prediction2(i,j,b)<DN_min
               prediction2_pos=max(DN_min,FRT1(i,j,b)+cubic_CR(i,j,b));
               prediction2(i,j,b)=min(DN_max,prediction2_pos);
            end
            if prediction2(i,j,b)>DN_max
               prediction2_pos=min(DN_max,FRT1(i,j,b)+cubic_CR(i,j,b));
               prediction2(i,j,b)=max(DN_min,prediction2_pos);
            end
            
        end
    end
end

%%%%%%%%%%%%%%%%%%%   Output    %%%%%%%%%%%%%%%%%%%%%%%%
disp('Outputting')

file_patch = prediction2;
info=enviinfo(file_patch);
enviwrite(file_patch,info,outputname);  

    
