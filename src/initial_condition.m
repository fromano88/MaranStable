function[b]=initial_condition(b, flowopt)
    energy=flowopt.energy;
    
    u_init=str2func(b.u_init.u);
    b.u.u=u_init(b.Z.u,b.R.u);
    b.u.w=interp2(b.ETA.u,b.XI.u,b.u.u,b.ETA.w,b.XI.w);
    b.u.p=interp2(b.ETA.u,b.XI.u,b.u.u,b.ETA.p,b.XI.p);
    b.u.T=interp2(b.ETA.u,b.XI.u,b.u.u,b.ETA.T,b.XI.T);
    b.u.v=interp2(b.ETA.u,b.XI.u,b.u.u,b.ETA.v,b.XI.v);
    
    w_init=str2func(b.w_init.w);
    b.w.w=w_init(b.Z.w,b.R_cyl.w);
    b.w.u=interp2(b.ETA.w,b.XI.w,b.w.w,b.ETA.u,b.XI.u);
    b.w.p=interp2(b.ETA.w,b.XI.w,b.w.w,b.ETA.p,b.XI.p);
    b.w.T=interp2(b.ETA.w,b.XI.w,b.w.w,b.ETA.T,b.XI.T);
    b.w.v=interp2(b.ETA.w,b.XI.w,b.w.w,b.ETA.v,b.XI.v);
    
    p_init=str2func(b.p_init.p);
    b.p.p=p_init(b.Z.p,b.R_cyl.p);
    b.p.v=interp1(b.xi.w(2:end-1).',b.p.p,b.xi.u.','linear','extrap');
    b.p.v=interp1(b.eta.u(2:end-1).',b.p.v.',b.eta.w.','linear','extrap').';
    
    T_init=str2func(b.T_init.T);
    b.T.T=T_init(b.Z.T,b.R_cyl.T);
    b.T.u=interp2(b.ETA.T,b.XI.T,b.T.T,b.ETA.u,b.XI.u);
    b.T.w=interp2(b.ETA.T,b.XI.T,b.T.T,b.ETA.w,b.XI.w);
    b.T.p=interp2(b.ETA.T,b.XI.T,b.T.T,b.ETA.p,b.XI.p);
    b.T.v=interp2(b.ETA.T,b.XI.T,b.T.T,b.ETA.v,b.XI.v);
    
    Rs1=[];
    Rsend=[];
    
    for n=1:length(b.bc.z)/2
        if strcmp(b.bc.z{n}(1),'s')==1 && strcmp(b.bc.z{n}(2),'s')==1 && eval(['b.' b.bc.z{n} '.o'])==1
            Rs1=[Rs1; b.R.u(1,(b.geom.z(n)<b.z.u & b.geom.z(n+1)>b.z.u)).'];
                            
        end
          
        if strcmp(b.bc.z{length(b.bc.z)/2+n}(1),'s')==1 && strcmp(b.bc.z{length(b.bc.z)/2+n}(2),'s')==1 && eval(['b.' b.bc.z{length(b.bc.z)/2+n} '.o'])==1
            Rsend=[Rsend; b.R.u(end,(b.geom.z(n)<b.z.u & b.geom.z(n+1)>b.z.u)).'];
                                    
        end
    end
    
    u=b.u.u.'; u=u(:);
    w=b.w.w.'; w=w(:);
    p=b.p.p.'; p=p(:);
    if energy==1
        T=b.T.T.'; T=T(:);
        b.x=[u;w;p;T;Rs1;Rsend];
        
    else
        b.x=[u;w;p;Rs1;Rsend];
        
    end