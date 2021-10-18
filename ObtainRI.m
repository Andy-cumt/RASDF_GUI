%%%%%%%%%%%%%%%%   Obtain RI  %%%%%%%%%%%%%%%%%
%CRT1:coarse image obtained at T1
%FRT1:fine image obtained at T1
%Factor: resolution ratio of coarse image to fine image

function RI =ObtainRI(CRT1,CRT2,FRT1,Factor)

%%%Resampling
[xH,yH,bands]=size(FRT1);
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
clear tmp1 tmp2 tmp3

nearest_FRT1C = (imresize(FRT1C,Factor,'nearest'));

%%%Linear regression
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

%%%Obtain RI
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
RI(RI<0.1)=0.1; %Avoid reliability index close to or equal to 0
