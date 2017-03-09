clear
close all
%%
load([pwd,'\simdata.mat']);
load([pwd,'\ExpParams.mat']);

%%
stack = stackread([exp.path,'\',exp.filename]);
stack = double(stack);
stack(find(stack(:)==0)) = 1;
[x,y,dxt,dyt,vxt,vyt] = cartPIV(stack,exp,pivspec,filtspec);
[dx,dy,vx,vy,vmag,ddiv,vdiv] = cartPIVstats(x,y,dxt,dyt,vxt,vyt);

%%
plotPIVoverlay(stack,exp,x,y,vx.t,vy.t,0,0,1,1,1,1.25,2.25)

