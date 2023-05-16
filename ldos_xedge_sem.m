(* ::Package:: *)

SetDirectory["~/Desktop/exciton"]
Nx = 6;
Ny = 32;
Ns = 500;
\[Eta]=0.004;
kernum=90;
mu = 0.1;
cdw = 0.111111;
wmax=0.25;
nw=500;
Nband=8;
flagx=1;
flagy=0;
ival=Table[i,{i,Nx*Ny*4/2}];
i=0;
For[iy=1,iy<=Ny,iy++,For[ix=1,ix<=Nx,ix++,i++;xpos[i]=ix;ypos[i]=iy;]];
NN=i;
gso[rxs_,rys_]:={{0.74DiracDelta[rxs]DiracDelta[rys]-(0.41+0.008I)DiracDelta[-1.+rxs,rys]-(0.41-0.008I)DiracDelta[1.+rxs,rys],0.14DiracDelta[1.-rxs,0.39+rys]-(0.+0.012I)DiracDelta[rxs,0.61-rys]-(0.+0.012I)DiracDelta[rxs,0.39+rys]-0.14DiracDelta[1.+rxs,0.39+rys],0.51DiracDelta[0.5-rxs,0.36-rys]+0.51DiracDelta[0.5+rxs,0.36-rys],0.39DiracDelta[0.5-rxs,0.25+rys]+0.29DiracDelta[1.5-rxs,0.25+rys]-0.39DiracDelta[0.5+rxs,0.25+rys]-0.29DiracDelta[1.5+rxs,0.25+rys],-0.031DiracDelta[1.-rxs,rys]+0.031DiracDelta[1.+rxs,rys],-0.05DiracDelta[rxs,0.61-rys]-0.051DiracDelta[rxs,0.39+rys],0,-0.011DiracDelta[0.5-rxs,0.25+rys]-0.011DiracDelta[0.5+rxs,0.25+rys]},{-0.14DiracDelta[1.-rxs,0.39-rys]+(0.+0.012I)DiracDelta[rxs,0.39-rys]+(0.+0.012I)DiracDelta[rxs,0.61+rys]+0.14DiracDelta[1.+rxs,0.39-rys],-1.6DiracDelta[rxs]DiracDelta[rys]+(1.13-0.01I)DiracDelta[-1.+rxs,rys]+0.13DiracDelta[rxs,1.-rys]+0.13DiracDelta[rxs,1.+rys]+(1.13+0.01I)DiracDelta[1.+rxs,rys],-0.39DiracDelta[0.5-rxs,0.25+rys]-0.29DiracDelta[1.5-rxs,0.25+rys]+0.39DiracDelta[0.5+rxs,0.25+rys]+0.29DiracDelta[1.5+rxs,0.25+rys],0.4DiracDelta[0.5-rxs,0.14-rys]+0.4DiracDelta[0.5+rxs,0.14-rys],0.051DiracDelta[rxs,0.39-rys]+0.05DiracDelta[rxs,0.61+rys],-0.04DiracDelta[1.-rxs,rys]+0.04DiracDelta[1.+rxs,rys],-0.011DiracDelta[0.5-rxs,0.25+rys]-0.011DiracDelta[0.5+rxs,0.25+rys],0},{0.51DiracDelta[0.5-rxs,0.36+rys]+0.51DiracDelta[0.5+rxs,0.36+rys],0.39DiracDelta[0.5-rxs,0.25-rys]+0.29DiracDelta[1.5-rxs,0.25-rys]-0.39DiracDelta[0.5+rxs,0.25-rys]-0.29DiracDelta[1.5+rxs,0.25-rys],0.74DiracDelta[rxs]DiracDelta[rys]-(0.41-0.008I)DiracDelta[-1.+rxs,rys]-(0.41+0.008I)DiracDelta[1.+rxs,rys],0.14DiracDelta[1.-rxs,0.39-rys]+(0.+0.012I)DiracDelta[rxs,0.39-rys]+(0.+0.012I)DiracDelta[rxs,0.61+rys]-0.14DiracDelta[1.+rxs,0.39-rys],0,0.011DiracDelta[0.5-rxs,0.25-rys]+0.011DiracDelta[0.5+rxs,0.25-rys],0.031DiracDelta[1.-rxs,rys]-0.031DiracDelta[1.+rxs,rys],0.051DiracDelta[rxs,0.39-rys]+0.05DiracDelta[rxs,0.61+rys]},{-0.39DiracDelta[0.5-rxs,0.25-rys]-0.29DiracDelta[1.5-rxs,0.25-rys]+0.39DiracDelta[0.5+rxs,0.25-rys]+0.29DiracDelta[1.5+rxs,0.25-rys],0.4DiracDelta[0.5-rxs,0.14+rys]+0.4DiracDelta[0.5+rxs,0.14+rys],-0.14DiracDelta[1.-rxs,0.39+rys]-(0.+0.012I)DiracDelta[rxs,0.61-rys]-(0.+0.012I)DiracDelta[rxs,0.39+rys]+0.14DiracDelta[1.+rxs,0.39+rys],-1.6DiracDelta[rxs]DiracDelta[rys]+(1.13+0.01I)DiracDelta[-1.+rxs,rys]+0.13DiracDelta[rxs,1.-rys]+0.13DiracDelta[rxs,1.+rys]+(1.13-0.01I)DiracDelta[1.+rxs,rys],0.011DiracDelta[0.5-rxs,0.25-rys]+0.011DiracDelta[0.5+rxs,0.25-rys],0,-0.05DiracDelta[rxs,0.61-rys]-0.051DiracDelta[rxs,0.39+rys],0.04DiracDelta[1.-rxs,rys]-0.04DiracDelta[1.+rxs,rys]},{0.031DiracDelta[1.-rxs,rys]-0.031DiracDelta[1.+rxs,rys],0.05DiracDelta[rxs,0.61-rys]+0.051DiracDelta[rxs,0.39+rys],0,0.011DiracDelta[0.5-rxs,0.25+rys]+0.011DiracDelta[0.5+rxs,0.25+rys],0.74DiracDelta[rxs]DiracDelta[rys]-(0.41-0.008I)DiracDelta[-1.+rxs,rys]-(0.41+0.008I)DiracDelta[1.+rxs,rys],0.14DiracDelta[1.-rxs,0.39+rys]+(0.+0.012I)DiracDelta[rxs,0.61-rys]+(0.+0.012I)DiracDelta[rxs,0.39+rys]-0.14DiracDelta[1.+rxs,0.39+rys],0.51DiracDelta[0.5-rxs,0.36-rys]+0.51DiracDelta[0.5+rxs,0.36-rys],0.39DiracDelta[0.5-rxs,0.25+rys]+0.29DiracDelta[1.5-rxs,0.25+rys]-0.39DiracDelta[0.5+rxs,0.25+rys]-0.29DiracDelta[1.5+rxs,0.25+rys]},{-0.051DiracDelta[rxs,0.39-rys]-0.05DiracDelta[rxs,0.61+rys],0.04DiracDelta[1.-rxs,rys]-0.04DiracDelta[1.+rxs,rys],0.011DiracDelta[0.5-rxs,0.25+rys]+0.011DiracDelta[0.5+rxs,0.25+rys],0,-0.14DiracDelta[1.-rxs,0.39-rys]-(0.+0.012I)DiracDelta[rxs,0.39-rys]-(0.+0.012I)DiracDelta[rxs,0.61+rys]+0.14DiracDelta[1.+rxs,0.39-rys],-1.6DiracDelta[rxs]DiracDelta[rys]+(1.13+0.01I)DiracDelta[-1.+rxs,rys]+0.13DiracDelta[rxs,1.-rys]+0.13DiracDelta[rxs,1.+rys]+(1.13-0.01I)DiracDelta[1.+rxs,rys],-0.39DiracDelta[0.5-rxs,0.25+rys]-0.29DiracDelta[1.5-rxs,0.25+rys]+0.39DiracDelta[0.5+rxs,0.25+rys]+0.29DiracDelta[1.5+rxs,0.25+rys],0.4DiracDelta[0.5-rxs,0.14-rys]+0.4DiracDelta[0.5+rxs,0.14-rys]},{0,-0.011DiracDelta[0.5-rxs,0.25-rys]-0.011DiracDelta[0.5+rxs,0.25-rys],-0.031DiracDelta[1.-rxs,rys]+0.031DiracDelta[1.+rxs,rys],-0.051DiracDelta[rxs,0.39-rys]-0.05DiracDelta[rxs,0.61+rys],0.51DiracDelta[0.5-rxs,0.36+rys]+0.51DiracDelta[0.5+rxs,0.36+rys],0.39DiracDelta[0.5-rxs,0.25-rys]+0.29DiracDelta[1.5-rxs,0.25-rys]-0.39DiracDelta[0.5+rxs,0.25-rys]-0.29DiracDelta[1.5+rxs,0.25-rys],0.74DiracDelta[rxs]DiracDelta[rys]-(0.41+0.008I)DiracDelta[-1.+rxs,rys]-(0.41-0.008I)DiracDelta[1.+rxs,rys],0.14DiracDelta[1.-rxs,0.39-rys]-(0.+0.012I)DiracDelta[rxs,0.39-rys]-(0.+0.012I)DiracDelta[rxs,0.61+rys]-0.14DiracDelta[1.+rxs,0.39-rys]},{-0.011DiracDelta[0.5-rxs,0.25-rys]-0.011DiracDelta[0.5+rxs,0.25-rys],0,0.05DiracDelta[rxs,0.61-rys]+0.051DiracDelta[rxs,0.39+rys],-0.04DiracDelta[1.-rxs,rys]+0.04DiracDelta[1.+rxs,rys],-0.39DiracDelta[0.5-rxs,0.25-rys]-0.29DiracDelta[1.5-rxs,0.25-rys]+0.39DiracDelta[0.5+rxs,0.25-rys]+0.29DiracDelta[1.5+rxs,0.25-rys],0.4DiracDelta[0.5-rxs,0.14+rys]+0.4DiracDelta[0.5+rxs,0.14+rys],-0.14DiracDelta[1.-rxs,0.39+rys]+(0.+0.012I)DiracDelta[rxs,0.61-rys]+(0.+0.012I)DiracDelta[rxs,0.39+rys]+0.14DiracDelta[1.+rxs,0.39+rys],-1.6DiracDelta[rxs]DiracDelta[rys]+(1.13-0.01I)DiracDelta[-1.+rxs,rys]+0.13DiracDelta[rxs,1.-rys]+0.13DiracDelta[rxs,1.+rys]+(1.13+0.01I)DiracDelta[1.+rxs,rys]}};
r[1]={0.,0.};
r[2]={0.,-0.39};
r[3]={0.5,-0.64};
r[4]={0.5,-0.25};
r[5]=r[1];r[6]=r[2];r[7]=r[3];r[8]=r[4];
For[m=1,m<=8,m++,For[n=1,n<=8,n++,For[cx=-2,cx<=2,cx++,For[cy=-2,cy<=2,cy++,rxs=cx+r[m][[1]]-r[n][[1]];rys=cy-r[m][[2]]+r[n][[2]];If[(cx==0&&cy==0&&m==n),t[-cx][-cy][m][n]=Chop[gso[rxs,rys][[m,n]]/DiracDelta[0]^2],t[-cx][-cy][m][n]=Chop[gso[rxs,rys][[m,n]]/DiracDelta[0,0]]];]]]];
Print["reading phasex"]
For[lx=1,lx<=NN*Nband,lx++,For[ly=1,ly<=NN*Nband,ly++,If[lx<=NN*Nband/2,m=Mod[lx-1/10,Nband/2]+1/10;i=Quotient[lx-0.1,Nband/2]+1;,m=Mod[lx-1/10,Nband/2]+1/10+4;i=Quotient[lx-0.1,Nband/2]+1-NN;];If[ly<=NN*Nband/2,n=Mod[ly-1/10,Nband/2]+1/10;j=Quotient[ly-0.1,Nband/2]+1;,n=Mod[ly-1/10,Nband/2]+1/10+4;j=Quotient[ly-0.1,Nband/2]+1-NN;];tijmn[lx][ly]=If[NumericQ[t[xpos[i]-(xpos[j]-Nx)][ypos[i]-(ypos[j]-Ny)][m][n]],-1.*flagx*flagy,If[NumericQ[t[xpos[i]-(xpos[j]-Nx)][ypos[i]-ypos[j]][m][n]],-1.*flagx,If[NumericQ[t[xpos[i]-(xpos[j]-Nx)][ypos[i]-(ypos[j]+Ny)][m][n]],-1.*flagx*flagy,If[NumericQ[t[xpos[i]-(xpos[j]+Nx)][ypos[i]-(ypos[j]-Ny)][m][n]],Plus[1.]*flagx*flagy,If[NumericQ[t[xpos[i]-(xpos[j]+Nx)][ypos[i]-ypos[j]][m][n]],Plus[1.]*flagx,If[NumericQ[t[xpos[i]-(xpos[j]+Nx)][ypos[i]-(ypos[j]+Ny)][m][n]],Plus[1.]*flagx*flagy]]]]]];]];
Phasex=SparseArray[{x_,y_}/;NumericQ[tijmn[x][y]]->tijmn[x][y],{NN*Nband,NN*Nband}];
Clear[tijmn]
Print["Reading H matrix now"]
h=Import["WTe2_Hamiltonian_Mtarix_CDW_sm_2pi6.wdx"]
Print["Reading H matrix done"]
Clear[fermi,spmtcalc]
Clear[spm,spmt,spmtab,spmij,SCspm,dijmn,dijmm]
Clear[kx,ky]
hsc=h*Exp[I*Re[kx]*Nx*Phasex];
Clear[hup,hud,hdu,hdn,rhup,ihup,rhud,ihud,rhdu,ihdu,rhdn,ihdn,h];
kkmx=Ns*Nx;
Print["starting ldos calc"]
ldos=0;
CloseKernels[];
LaunchKernels[kernum];
ldos=Total[ParallelTable[ES=Eigensystem[hsc];elist=Table[ES[[1,l]],{l,1,Nband*NN}];ulist=Table[Abs[ES[[2,l,ival]]]^2,{l,1,Nband*NN}];ldossc=Table[Im[Total[Table[ulist[[l,1;;All]]*(1/(-wmax+wmax*2*(w/nw)-elist[[l]]+I*\[Eta])),{l,1,Nband*NN}]]],{w,0,nw}];ldossc,{kx,0,(Ns-1.)*2*(Pi/kkmx),2*(Pi/kkmx)}]];
Print["done"];
fl=StringJoin["xedge_sem_cdw",ToString[cdw],"_mu",ToString[mu],"_wmax",ToString[wmax],"Ny",ToString[Ny],"Nx",ToString[Nx],"eta",ToString[\[Eta]],"Ns",ToString[Ns]];
DeleteFile[StringJoin[fl,"Ldos_wte2.dat"]];
str1=OpenAppend[StringJoin[fl,"Ldos_wte2.dat"],FormatType->StandardForm];
For[n=1,n<=nw+1,n++,For[i=1,i<=Length[ival],i++,Write[str1,-wmax+wmax*2*((n-1)/nw),"  ",(-1/(Pi))*ldos[[n,i]]/(Ns)];]];
Close[str1];
Exit[];
