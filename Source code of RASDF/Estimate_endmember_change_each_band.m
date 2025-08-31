
function endmember_change = Estimate_endmember_change_each_band(fraction_T1,CR_T2,CR_T1,xL,yL,class,bands,mm,nn,mmm,nnn,min_allowb,max_allowb,p,pp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_change=length(mm(:));
if p>=2
    for ppp=1:num_change
        fraction_T1(mm(ppp),nn(ppp),:)=0;
    end
end
if pp>=2
    for ppp=1:pp-1
        fraction_T1(mmm(ppp),nnn(ppp),:)=0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%find the number of candidate pixel for unmixing
number=0;
for ii=1:xL
    for jj=1:yL
        if sum(fraction_T1(ii,jj,:))~=0
            number=number+1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if number<=class*4
    endmember_change = zeros(bands,class);
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    location_fraction=zeros(xL,yL,class);
    location=zeros(xL,yL,class);
    k=100;



    endmember_change = zeros(bands,class);
    np_estimate = zeros(bands,1);

    for b=1:bands

        y_delta = CR_T2(:,:,b)-CR_T1(:,:,b);
        min_allow=min_allowb - std(reshape(y_delta,xL*yL,1));
        max_allow=max_allowb + std(reshape(y_delta,xL*yL,1));
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for c=1:class

        A=reshape(fraction_T1(:,:,c),xL*yL,1);
        [fraction_sort, D]= sort(A,'descend');
        for n=1:k
            x_tmp=mod(D(n),xL);
            y_tmp=ceil(D(n)/xL);
            if x_tmp==0
                x_tmp=xL;
            end
            if fraction_sort(n)>0.01 % 
                location_fraction(x_tmp,y_tmp,c) = 1;
            end
        end

        num_f = sum(sum(fraction_T1(:,:,c)>0.01));    
        num_pure1 = min(k,num_f);    
        y_tmp = abs(y_delta).*(location_fraction(:,:,c)==1) + (location_fraction(:,:,c)~=1)*999999999999999999;
        A=reshape(y_tmp,xL*yL,1);
        [y_sort, D]= sort(A,'ascend');
        tmp = y_sort(1:num_pure1);
        tmp1 = tmp(ceil(0.01*num_pure1));
        tmp2 = tmp(ceil(1.0*num_pure1));
        location(:,:,c) = (location_fraction(:,:,c)==1) .* double(y_tmp>=tmp1) .* double(y_tmp<=tmp2) ;
    end

    thresholdMatrix=zeros(xL,yL);
    for c=1:class
        thresholdMatrix = thresholdMatrix+location(:,:,c);
    end
    thresholdMatrix=thresholdMatrix>0;

    np_estimate(b) = sum(sum(thresholdMatrix));

        x0 = ones(1,class) / class;    
        lb = ones(1,class) * min_allow;
        ub = ones(1,class) * max_allow;            
        options = optimoptions('fmincon','Display','off');

        [x,fval] = fmincon(@(x) objfun_deltaEndm(x,fraction_T1,y_delta,thresholdMatrix,xL,yL),x0,[],[],[],[],lb,ub,[],options);
        endmember_change(b,:) = x;

    end
end
 

