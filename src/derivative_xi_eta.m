function drv = derivative_xi_eta(b,var,dir,loc)


switch dir
    case 'xi'
        if strcmp(loc,'p')
            if isfield(var,'u')
                drv = (var.u(2:end,2:end-1)-var.u(1:end-1,2:end-1))./b.DXI.c;
            elseif isfield(var,'w')
                drv = (var.v(2:end,:)-var.v(1:end-1,:))./repmat(b.dxi.c',1,size(b.R.v,2));
                drv = interp2(b.ETA.w(2:end-1,:),b.XI.w(2:end-1,:),drv,b.ETA.p,b.XI.p,'spline');
            else
                drv = (var.T(2:end,2:end-1)-var.T(1:end-1,2:end-1))./b.DXI.rm;
                drv = interp2(b.ETA.u(:,2:end-1),b.XI.u(:,2:end-1),drv,b.ETA.p,b.XI.p,'spline');
            end
            
        elseif strcmp(loc,'u')
            drv = (var.T(2:end,2:end-1)-var.T(1:end-1,2:end-1))./b.DXI.rm;
            
        elseif strcmp(loc,'w')
            drv = (var.v(2:end,:)-var.v(1:end-1,:))./repmat(b.dxi.c',1,size(b.R.v,2));
            
        elseif strcmp(loc,'v')
            drv = (var.w(2:end,:)-var.w(1:end-1,:))./repmat(b.dxi.rm',1,size(b.R.v,2));
            
        end
        
        
    case 'eta'
        if strcmp(loc,'p')
            if isfield(var,'w')
                drv = (var.w(2:end-1,2:end)-var.w(2:end-1,1:end-1))./b.DETA.c;
            elseif isfield(var,'u')
                drv = (var.v(:,2:end)-var.v(:,1:end-1))./repmat(b.deta.c,size(b.R.v,1),1);
                drv = interp2(b.ETA.u(:,2:end-1),b.XI.u(:,2:end-1),drv,b.ETA.p,b.XI.p,'spline');
            else
                drv = (var.T(2:end-1,2:end)-var.T(2:end-1,1:end-1))./b.DETA.zm;
                drv = interp2(b.ETA.w(2:end-1,:),b.XI.w(2:end-1,:),drv,b.ETA.p,b.XI.p,'spline');
            end
            
        elseif strcmp(loc,'w')
            drv = (var.T(2:end-1,2:end)-var.T(2:end-1,1:end-1))./b.DETA.zm;
            
        elseif strcmp(loc,'v')
            drv = (var.u(:,2:end)-var.u(:,1:end-1))./repmat(b.deta.zm,size(b.R.v,1),1);
            
        elseif strcmp(loc,'h') % on the surface
            drv = (var.w(end,2:end)-var.w(end,1:end-1))./b.DETA.c(end,:);
            
        end
            
        
    otherwise
        error('Choose either derivative in ''xi'' or ''eta'' direction!')
        
end