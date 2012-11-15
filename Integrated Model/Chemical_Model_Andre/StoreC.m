function [status] = StoreC(t,y,x)

global C C2 Vv Call Call2 V2 tt  
if strcmp(x,'done')==0 && strcmp(x,'init')==0
     
    Call  = [Call C];
    Call2 = [Call2 C2];
    V2 = [V2 Vv];
    tt = [tt;t(end)];
end
status = 0;
end

