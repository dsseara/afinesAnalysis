function [dx,dy,vx,vy,vmag,ddiv,vdiv] = cartPIVstats(x,y,dxt,dyt,vxt,vyt)
dim = size(dxt,3);

dx.t = dxt;
dx.mean = nanmean(reshape(dxt,[size(dxt,1)*size(dxt,2), size(dxt,3)]));
dx.std = nanstd(reshape(dxt,[size(dxt,1)*size(dxt,2), size(dxt,3)]));

dy.t = dyt;
dy.mean = nanmean(reshape(dyt,[size(dyt,1)*size(dyt,2), size(dyt,3)]));
dy.std = nanstd(reshape(dyt,[size(dyt,1)*size(dyt,2), size(dyt,3)]));

vx.t = vxt;
vx.mean = nanmean(reshape(vxt,[size(vxt,1)*size(vxt,2), size(vxt,3)]));
vx.std = nanstd(reshape(vxt,[size(vxt,1)*size(vxt,2), size(vxt,3)]));

vy.t = vyt;
vy.mean = nanmean(reshape(vyt,[size(vyt,1)*size(vyt,2), size(vyt,3)]));
vy.std = nanstd(reshape(vyt,[size(vyt,1)*size(vyt,2), size(vyt,3)]));

vmagt = sqrt(vxt.^2+vyt.^2);
vmag.t = vmagt;
vmag.mean = nanmean(reshape(vmagt,[size(vmagt,1)*size(vmagt,2), size(vmagt,3)]));
vmag.std = nanstd(reshape(vmagt,[size(vmagt,1)*size(vmagt,2), size(vmagt,3)]));

ddivt = NaN(size(dxt));
vdivt = NaN(size(vxt));
for i=1:dim
    ddivt(:,:,i) = divergence(x,y,dxt(:,:,i),dyt(:,:,i));
    vdivt(:,:,i) = divergence(x,y,vxt(:,:,i),vyt(:,:,i));
end

ddiv.t = ddivt;
ddiv.mean = nanmean(reshape(ddivt,[size(ddivt,1)*size(ddivt,2), size(ddivt,3)]));
ddiv.std = nanstd(reshape(ddivt,[size(ddivt,1)*size(ddivt,2), size(ddivt,3)]));

vdiv.t = vdivt;
vdiv.mean = nanmean(reshape(vdivt,[size(vdivt,1)*size(vdivt,2), size(vdivt,3)]));
vdiv.std = nanstd(reshape(vdivt,[size(vdivt,1)*size(vdivt,2), size(vdivt,3)]));

clearvars -except dx dy vx vy vmag ddiv vdiv

if exist([pwd,'/WS'], 'dir') == 0, mkdir([pwd,'/WS']); end
save([pwd,'/WS/cartPIVstats_WS.mat']);
end

