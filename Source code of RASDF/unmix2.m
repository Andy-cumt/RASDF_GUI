function endmember_change = unmix2(fraction_T1,CR_T2,CR_T1,xL,yL,class,bands)

endmember_change = zeros(bands,class);

for b=1:bands
    y_delta = (CR_T2(:,b)-CR_T1(:,b));
    min_allow = min(reshape(y_delta,xL*yL,1)) - std(reshape(y_delta,xL*yL,1));
    max_allow = max(reshape(y_delta,xL*yL,1)) + std(reshape(y_delta,xL*yL,1));
    x0 = ones(1,class) / class;    
    lb = ones(1,class) * min_allow;
    ub = ones(1,class) * max_allow;            
    options = optimoptions('fmincon','Display','off');
    thresholdMatrix=ones(xL,yL); 
    [x,fval] = fmincon(@(x) objfun_deltaEndm(x,fraction_T1,y_delta,thresholdMatrix,xL,yL),x0,[],[],[],[],lb,ub,[],options);
    endmember_change(b,:) = x;
end
 

