function [xn,err,xn_b,err_b,err_x,err_x_b] = defiltering_Steffensen(ys, f, gammab, early_stop, max_iter, xs)
if size(ys, 3) == 1
    [xn,err,xn_b,err_b,err_x,err_x_b] = defiltering_Steffensen1(ys, f, gammab, early_stop, max_iter,xs);
else
    r = ys(:,:,1);
    g = ys(:,:,2);
    b = ys(:,:,3);
    
    if (size(xs,1)==1)
        xr = 0;
        xg = 0;
        xb = 0;
    else
        xr = xs(:,:,1);
        xg = xs(:,:,2);
        xb = xs(:,:,3);
    end
    
    [xnr,errr,xnr_b,errr_b,errr_x,errr_x_b] = defiltering_Steffensen1(r, f, gammab, early_stop, max_iter, xr);
    [xng,errg,xng_b,errg_b,errg_x,errg_x_b] = defiltering_Steffensen1(g, f, gammab, early_stop, max_iter, xg);
    [xnb,errb,xnb_b,errb_b,errb_x,errb_x_b] = defiltering_Steffensen1(b, f, gammab, early_stop, max_iter, xb);
    
    xn = cat(3, xnr, xng, xnb);
    xn_b = cat(3, xnr_b, xng_b, xnb_b);
    
    if (early_stop==0)
        err = (errr + errg + errb)./3;
        err_b = (errr_b + errg_b + errb_b)./3;
        err_x = (errr_x + errg_x + errb_x)./3;
        err_x_b = (errr_x_b + errg_x_b + errb_x_b)./3;
    else
        % take the longest error vector
        lr = length(errr);
        lg = length(errg);
        lb = length(errb);
        m = max([lr, lg, lb]);
        if m==lr, err=errr; err_b=errr_b; err_x=errr_x; err_x_b=errr_x_b; end
        if m==lg, err=errg; err_b=errg_b; err_x=errg_x; err_x_b=errg_x_b; end
        if m==lb, err=errb; err_b=errb_b; err_x=errb_x; err_x_b=errb_x_b; end
    end
    
end
end

function [xn,err,xn_b,err_b,err_x,err_x_b] = defiltering_Steffensen1(ys, f, gammab, early_stop, max_iter, xs)
% Steffensen xn=xn+f(xn)/g(xn), g(x)=(f(x+h)-f(x))/h, h=y-f(x)

xn = 2*ys-f(ys);
%xn = ys;

err = [];
fxn = f(xn);
ysfxn = ys-fxn;
nys = norm(ys(:));
err = [err; norm(ysfxn(:))/nys];

err_b = [];
err_b = [err_b; norm(ysfxn(:))/nys];
xn_b = xn;

err_x = [];
if (size(xs,1)~=1)
    xsxn = xs-xn;
    nxs=norm(xs(:));
    err_x = [err_x; norm(xsxn(:))/nxs];
end
err_x_b = [];
if (size(xs,1)~=1)
    xsxn = xs-xn;
    nxs=norm(xs(:));
    err_x_b = [err_x_b; norm(xsxn(:))/nxs];
end

for i=1:max_iter
    xn = xn + norm(ysfxn)*ysfxn/norm(f(xn+ysfxn)-f(xn));
        
    fxn = f(xn);
    ysfxn = ys-fxn;
    
    errn = norm(ysfxn(:))/nys;
    err = [err; errn];
    errn_b = err_b(length(err_b));
    err_b = [err_b; min(errn_b, errn)];
    if (errn < errn_b)
        xn_b = xn;
    end
    
    if (size(xs,1)~=1)
        xsxn = xs-xn;
        err_x = [err_x; norm(xsxn(:))/nxs];
        xsxnb = xs-xn_b;
        err_x_b = [err_x_b; norm(xsxnb(:))/nxs];
    end
    
    if (err(length(err))<1e-3) && (early_stop==1)
        fprintf('Break at iter: %d || Error: %f\n', i, err(length(err)));
        break;
    end
end
end
