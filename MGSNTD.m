function [G,A,B,C,iter]=MGSNTD2024(X,G,A,B,C,LW,PF,PY,c,options)
% NTD+lambda*Tr(C'LC)+lambda*Tr(C'L~C)+beta*||PF-PY||_(F)^(2)
% LW: the sum of weighting matrices W and W~
% c: orders of labeled samples 

eps=options.eps;
Maxiter=options.Maxiter;
lambda=options.lambda;
beta=options.beta;
LD=sum(LW,2);LD=diag(LD);   
PY(c,:)=C(c,:);
PY= sparse(PY);
for k=1:Maxiter
    %% 更新U  
    G1=tenmat(G,1);
    G1=G1.data;
    P=ttm(X,{B',C'},[2,3]);
    P1=tenmat(P,1);
    P1=P1.data;
    AFZ=P1*G1';
    %YWVS=max(YWVS,0);  
    AFZ(AFZ<0)=0;
    QQ=ttm(G,{B'*B,C'*C},[2,3]);
    Q1=tenmat(QQ,1);
    Q1=Q1.data;
    AFM=A*Q1*G1';  
    A=A.*(AFZ./max(AFM,1e-15));
    %U=UU;
    %% 更新V
    G2=tenmat(G,2);
    G2=G2.data;
    R=ttm(X,{A',C'},[1,3]);
    R2=tenmat(R,2);
    R2=R2.data;
    BFZ=R2*G2';
    %YWUS=max(YWUS,0); 
    BFZ(BFZ<0)=0;
    E=ttm(G,{A'*A,C'*C},[1,3]);
    E2=tenmat(E,2);
    E2=E2.data;
    BFM=B*E2*G2';   
    B=B.*(BFZ./max(BFM,1e-15));
    %V=VV;
    %% 更新W
    W0=C;
    G3=tenmat(G,3);
    G3=G3.data;
    F=ttm(X,{A',B'},[1,2]);
    F3=tenmat(F,3);
    F3=F3.data;
    YVUS=F3*G3';
    CFZ=YVUS+lambda*LW*C+beta*PF;%LS=W(S),LD=D,L=D-W(S)
    CFZ(CFZ<0)=0;
   %YL=max(YL,0);
    H=ttm(G,{A'*A,B'*B},[1,2]);
    H3=tenmat(H,3);
    H3=H3.data;
    CFM1=C*H3*G3';
    la=lambda*LD*C+beta*PY;
    CFM=CFM1+la;
    C=C.*(CFZ./max(CFM,1e-15));
    PY(c,:)=C(c,:);
    %% 更新G
    GFZ=ttm(X,{A',B',C'},[1,2,3]);
    GFZ=GFZ.data;
    GFM=ttm(G,{A'*A,B'*B,C'*C},[1,2,3]);
    GFM=GFM.data;
    G4=G.data;
    G4=G4.*(GFZ./max(GFM,1e-15));
    G=tensor(G4);
    %S=tensor(reshape(S_vec, dimU, dimV, dimW));   
   %%停止准则
    Re=C-W0;
    Err=norm(Re(:))/norm(C(:));
    iter=k;
    if Err<eps
        break        
    end     
 
end
end