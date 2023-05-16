(* ::Package:: *)

SetDirectory["~/Desktop/exciton/4_orb"];


Nx=80;
Ny=6;
Ns=4000;
\[Eta] = 0.002;
kernum=42;
wmax = 0.25;
mu = -0.25;
cdw = -0.08;
nw = 500; 
Nband=8;
lyr = Nx/2;
flagx=0;
flagy=1;
ival=Flatten[Table[{j+i,j+i+1,j+i+2,j+i+3},{j,0,(lyr-1)*4,Nband/2},{i,1,Ny*Nx*Nband/2,Nx*Nband/2}]];
i=0;
For[iy = 1, iy <= Ny, iy++,
For[ix = 1, ix <= Nx, ix++, i++;
    xpos[i] = ix;
    ypos[i] = iy;
    ]];
      NN=i;


gso[rxs_,rys_]:={{0.74 DiracDelta[rxs] DiracDelta[rys]-(0.41+0.008 I) DiracDelta[-1.+rxs,rys]-(0.41-0.008 I) DiracDelta[1.+rxs,rys],0.14 DiracDelta[1.-rxs,0.39+rys]-(0.+0.012 I) DiracDelta[rxs,0.61-rys]-(0.+0.012 I) DiracDelta[rxs,0.39+rys]-0.14 DiracDelta[1.+rxs,0.39+rys],0.51 DiracDelta[0.5-rxs,0.36-rys]+0.51 DiracDelta[0.5+rxs,0.36-rys],0.39 DiracDelta[0.5-rxs,0.25+rys]+0.29 DiracDelta[1.5-rxs,0.25+rys]-0.39 DiracDelta[0.5+rxs,0.25+rys]-0.29 DiracDelta[1.5+rxs,0.25+rys],-0.031 DiracDelta[1.-rxs,rys]+0.031 DiracDelta[1.+rxs,rys],-0.05 DiracDelta[rxs,0.61-rys]-0.051 DiracDelta[rxs,0.39+rys],0,-0.011 DiracDelta[0.5-rxs,0.25+rys]-0.011 DiracDelta[0.5+rxs,0.25+rys]},{-0.14 DiracDelta[1.-rxs,0.39-rys]+(0.+0.012 I) DiracDelta[rxs,0.39-rys]+(0.+0.012 I) DiracDelta[rxs,0.61+rys]+0.14 DiracDelta[1.+rxs,0.39-rys],-1.75 DiracDelta[rxs] DiracDelta[rys]+(1.13-0.01 I) DiracDelta[-1.+rxs,rys]+0.13 DiracDelta[rxs,1.-rys]+0.13 DiracDelta[rxs,1.+rys]+(1.13+0.01 I) DiracDelta[1.+rxs,rys],-0.39 DiracDelta[0.5-rxs,0.25+rys]-0.29 DiracDelta[1.5-rxs,0.25+rys]+0.39 DiracDelta[0.5+rxs,0.25+rys]+0.29 DiracDelta[1.5+rxs,0.25+rys],0.4 DiracDelta[0.5-rxs,0.14-rys]+0.4 DiracDelta[0.5+rxs,0.14-rys],0.051 DiracDelta[rxs,0.39-rys]+0.05 DiracDelta[rxs,0.61+rys],-0.04 DiracDelta[1.-rxs,rys]+0.04 DiracDelta[1.+rxs,rys],-0.011 DiracDelta[0.5-rxs,0.25+rys]-0.011 DiracDelta[0.5+rxs,0.25+rys],0},{0.51 DiracDelta[0.5-rxs,0.36+rys]+0.51 DiracDelta[0.5+rxs,0.36+rys],0.39 DiracDelta[0.5-rxs,0.25-rys]+0.29 DiracDelta[1.5-rxs,0.25-rys]-0.39 DiracDelta[0.5+rxs,0.25-rys]-0.29 DiracDelta[1.5+rxs,0.25-rys],0.74 DiracDelta[rxs] DiracDelta[rys]-(0.41-0.008 I) DiracDelta[-1.+rxs,rys]-(0.41+0.008 I) DiracDelta[1.+rxs,rys],0.14 DiracDelta[1.-rxs,0.39-rys]+(0.+0.012 I) DiracDelta[rxs,0.39-rys]+(0.+0.012 I) DiracDelta[rxs,0.61+rys]-0.14 DiracDelta[1.+rxs,0.39-rys],0,0.011 DiracDelta[0.5-rxs,0.25-rys]+0.011 DiracDelta[0.5+rxs,0.25-rys],0.031 DiracDelta[1.-rxs,rys]-0.031 DiracDelta[1.+rxs,rys],0.051 DiracDelta[rxs,0.39-rys]+0.05 DiracDelta[rxs,0.61+rys]},{-0.39 DiracDelta[0.5-rxs,0.25-rys]-0.29 DiracDelta[1.5-rxs,0.25-rys]+0.39 DiracDelta[0.5+rxs,0.25-rys]+0.29 DiracDelta[1.5+rxs,0.25-rys],0.4 DiracDelta[0.5-rxs,0.14+rys]+0.4 DiracDelta[0.5+rxs,0.14+rys],-0.14 DiracDelta[1.-rxs,0.39+rys]-(0.+0.012 I) DiracDelta[rxs,0.61-rys]-(0.+0.012 I) DiracDelta[rxs,0.39+rys]+0.14 DiracDelta[1.+rxs,0.39+rys],-1.75 DiracDelta[rxs] DiracDelta[rys]+(1.13+0.01 I) DiracDelta[-1.+rxs,rys]+0.13 DiracDelta[rxs,1.-rys]+0.13 DiracDelta[rxs,1.+rys]+(1.13-0.01 I) DiracDelta[1.+rxs,rys],0.011 DiracDelta[0.5-rxs,0.25-rys]+0.011 DiracDelta[0.5+rxs,0.25-rys],0,-0.05 DiracDelta[rxs,0.61-rys]-0.051 DiracDelta[rxs,0.39+rys],0.04 DiracDelta[1.-rxs,rys]-0.04 DiracDelta[1.+rxs,rys]},{0.031 DiracDelta[1.-rxs,rys]-0.031 DiracDelta[1.+rxs,rys],0.05 DiracDelta[rxs,0.61-rys]+0.051 DiracDelta[rxs,0.39+rys],0,0.011 DiracDelta[0.5-rxs,0.25+rys]+0.011 DiracDelta[0.5+rxs,0.25+rys],0.74 DiracDelta[rxs] DiracDelta[rys]-(0.41-0.008 I) DiracDelta[-1.+rxs,rys]-(0.41+0.008 I) DiracDelta[1.+rxs,rys],0.14 DiracDelta[1.-rxs,0.39+rys]+(0.+0.012 I) DiracDelta[rxs,0.61-rys]+(0.+0.012 I) DiracDelta[rxs,0.39+rys]-0.14 DiracDelta[1.+rxs,0.39+rys],0.51 DiracDelta[0.5-rxs,0.36-rys]+0.51 DiracDelta[0.5+rxs,0.36-rys],0.39 DiracDelta[0.5-rxs,0.25+rys]+0.29 DiracDelta[1.5-rxs,0.25+rys]-0.39 DiracDelta[0.5+rxs,0.25+rys]-0.29 DiracDelta[1.5+rxs,0.25+rys]},{-0.051 DiracDelta[rxs,0.39-rys]-0.05 DiracDelta[rxs,0.61+rys],0.04 DiracDelta[1.-rxs,rys]-0.04 DiracDelta[1.+rxs,rys],0.011 DiracDelta[0.5-rxs,0.25+rys]+0.011 DiracDelta[0.5+rxs,0.25+rys],0,-0.14 DiracDelta[1.-rxs,0.39-rys]-(0.+0.012 I) DiracDelta[rxs,0.39-rys]-(0.+0.012 I) DiracDelta[rxs,0.61+rys]+0.14 DiracDelta[1.+rxs,0.39-rys],-1.75 DiracDelta[rxs] DiracDelta[rys]+(1.13+0.01 I) DiracDelta[-1.+rxs,rys]+0.13 DiracDelta[rxs,1.-rys]+0.13 DiracDelta[rxs,1.+rys]+(1.13-0.01 I) DiracDelta[1.+rxs,rys],-0.39 DiracDelta[0.5-rxs,0.25+rys]-0.29 DiracDelta[1.5-rxs,0.25+rys]+0.39 DiracDelta[0.5+rxs,0.25+rys]+0.29 DiracDelta[1.5+rxs,0.25+rys],0.4 DiracDelta[0.5-rxs,0.14-rys]+0.4 DiracDelta[0.5+rxs,0.14-rys]},{0,-0.011 DiracDelta[0.5-rxs,0.25-rys]-0.011 DiracDelta[0.5+rxs,0.25-rys],-0.031 DiracDelta[1.-rxs,rys]+0.031 DiracDelta[1.+rxs,rys],-0.051 DiracDelta[rxs,0.39-rys]-0.05 DiracDelta[rxs,0.61+rys],0.51 DiracDelta[0.5-rxs,0.36+rys]+0.51 DiracDelta[0.5+rxs,0.36+rys],0.39 DiracDelta[0.5-rxs,0.25-rys]+0.29 DiracDelta[1.5-rxs,0.25-rys]-0.39 DiracDelta[0.5+rxs,0.25-rys]-0.29 DiracDelta[1.5+rxs,0.25-rys],0.74 DiracDelta[rxs] DiracDelta[rys]-(0.41+0.008 I) DiracDelta[-1.+rxs,rys]-(0.41-0.008 I) DiracDelta[1.+rxs,rys],0.14 DiracDelta[1.-rxs,0.39-rys]-(0.+0.012 I) DiracDelta[rxs,0.39-rys]-(0.+0.012 I) DiracDelta[rxs,0.61+rys]-0.14 DiracDelta[1.+rxs,0.39-rys]},{-0.011 DiracDelta[0.5-rxs,0.25-rys]-0.011 DiracDelta[0.5+rxs,0.25-rys],0,0.05 DiracDelta[rxs,0.61-rys]+0.051 DiracDelta[rxs,0.39+rys],-0.04 DiracDelta[1.-rxs,rys]+0.04 DiracDelta[1.+rxs,rys],-0.39 DiracDelta[0.5-rxs,0.25-rys]-0.29 DiracDelta[1.5-rxs,0.25-rys]+0.39 DiracDelta[0.5+rxs,0.25-rys]+0.29 DiracDelta[1.5+rxs,0.25-rys],0.4 DiracDelta[0.5-rxs,0.14+rys]+0.4 DiracDelta[0.5+rxs,0.14+rys],-0.14 DiracDelta[1.-rxs,0.39+rys]+(0.+0.012 I) DiracDelta[rxs,0.61-rys]+(0.+0.012 I) DiracDelta[rxs,0.39+rys]+0.14 DiracDelta[1.+rxs,0.39+rys],-1.75 DiracDelta[rxs] DiracDelta[rys]+(1.13-0.01 I) DiracDelta[-1.+rxs,rys]+0.13 DiracDelta[rxs,1.-rys]+0.13 DiracDelta[rxs,1.+rys]+(1.13+0.01 I) DiracDelta[1.+rxs,rys]}};


r[1]={0.,0.};
r[2]={0.,-0.39};
r[3]={0.5,-0.64};
r[4]={0.5,-0.25};
r[5]=r[1];r[6]=r[2];r[7]=r[3];r[8]=r[4];
For[m=1,m<=8,m++,
For[n=1,n<=8,n++,
For[cx=-2,cx<=2,cx++,
For[cy=-2,cy<=2,cy++,
rxs=cx+r[m][[1]]-r[n][[1]];rys=cy-r[m][[2]]+r[n][[2]];
If[(cx==0 && cy==0 && m==n),t[-cx][-cy][m][n]=Chop[gso[rxs,rys][[m,n]]/DiracDelta[0]^2],
t[-cx][-cy][m][n]=Chop[gso[rxs,rys][[m,n]]/DiracDelta[0,0]]];
]]]];


(*Print["reading phasex"]*)


(*For[lx=1,lx<=NN*Nband,lx++,
For[ly=1,ly<=NN*Nband,ly++,
If[lx<= NN*Nband/2,
m=Mod[lx-1/10,Nband/2]+1/10;i=Quotient[lx-0.1,Nband/2]+1;,
m=Mod[lx-1/10,Nband/2]+1/10+4;i=Quotient[lx-0.1,Nband/2]+1-NN;];
If[ly<= NN*Nband/2,
n=Mod[ly-1/10,Nband/2]+1/10;j=Quotient[ly-0.1,Nband/2]+1;,
n=Mod[ly-1/10,Nband/2]+1/10+4;j=Quotient[ly-0.1,Nband/2]+1-NN;];
tijmn[lx][ly]=If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] - Ny)][m][n]], -1.*flagx*flagy, 
If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - ypos[j]][m][n]], -1.*flagx, 
If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] + Ny)][m][n]], -1.*flagx*flagy, 
If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] - Ny)][m][n]],Plus[1.]*flagx*flagy, 
If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - ypos[j]][m][n]], Plus[1.]*flagx, 
If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] + Ny)][m][n]], Plus[1.]*flagx*flagy]]]]]];
]];
Phasex=SparseArray[{x_,y_}/;NumericQ[tijmn[x][y]]->tijmn[x][y],{NN*Nband,NN*Nband}];
Clear[tijmn]*)


(*Phasexup=Take[Phasex,{1,NN*Nband/2},{1,NN*Nband/2}];
Phasexud=Take[Phasex,{1,NN*Nband/2},{NN*Nband/2+1,NN*Nband}];
Phasexdu=Take[Phasex,{NN*Nband/2+1,NN*Nband},{1,NN*Nband/2}];
Phasexdn=Take[Phasex,{NN*Nband/2+1,NN*Nband},{NN*Nband/2+1,NN*Nband}];
Clear[Phasex];*)


Print["reading phasey"]


For[lx=1,lx<=NN*Nband,lx++,
For[ly=1,ly<=NN*Nband,ly++,
If[lx<= NN*Nband/2,
m=Mod[lx-1/10,Nband/2]+1/10;i=Quotient[lx-0.1,Nband/2]+1;,
m=Mod[lx-1/10,Nband/2]+1/10+4;i=Quotient[lx-0.1,Nband/2]+1-NN;];
If[ly<= NN*Nband/2,
n=Mod[ly-1/10,Nband/2]+1/10;j=Quotient[ly-0.1,Nband/2]+1;,
n=Mod[ly-1/10,Nband/2]+1/10+4;j=Quotient[ly-0.1,Nband/2]+1-NN;];
tijmn[lx][ly]=If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] - Ny)][m][n]], -1.*flagx*flagy,
 If[NumericQ[t[xpos[i] - (xpos[j] - Nx)][ypos[i] - (ypos[j] + Ny)][m][n]], Plus[1.]*flagx*flagy,
  If[NumericQ[t[xpos[i] - xpos[j]][ypos[i] - (ypos[j] - Ny)][m][n]], -1.*flagy, 
        If[NumericQ[t[xpos[i] - xpos[j]][ypos[i] - (ypos[j] + Ny)][m][n]], Plus[1.]*flagy,
         If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] - Ny)][m][n]], -1.*flagx*flagy,
          If[NumericQ[t[xpos[i] - (xpos[j] + Nx)][ypos[i] - (ypos[j] + Ny)][m][n]], Plus[1.]*flagx*flagy]]]]]];
          ]]; 
Phasey=SparseArray[{x_,y_}/;NumericQ[tijmn[x][y]]->tijmn[x][y],{NN*Nband,NN*Nband}];
Clear[tijmn]


(*Phaseyup=Take[Phasey,{1,NN*Nband/2},{1,NN*Nband/2}];
Phaseyud=Take[Phasey,{1,NN*Nband/2},{NN*Nband/2+1,NN*Nband}];
Phaseydu=Take[Phasey,{NN*Nband/2+1,NN*Nband},{1,NN*Nband/2}];
Phaseydn=Take[Phasey,{NN*Nband/2+1,NN*Nband},{NN*Nband/2+1,NN*Nband}];
Clear[Phasey];*)


Print["Reading H matrix now"]
H = Import["WTe2_Hamiltonian_Mtarix_CDW_yedge_sem_2pi6.wdx"]
(*rhup = Import["rhup.dat"];
ihup = Import["ihup.dat"];
rhud = Import["rhud.dat"];
ihud = Import["ihud.dat"];
rhdu = Import["rhdu.dat"];
ihdu = Import["ihdu.dat"];
rhdn = Import["rhdn.dat"];
ihdn = Import["ihdn.dat"];*)

Print["Reading H matrix done"]


(*hup=rhup+I*ihup;
hud=rhud+I*ihud;
hdu=rhdu+I*ihdu;
hdn=rhdn+I*ihdn;
*)


Clear[fermi, spmtcalc]
Clear[spm, spmt, spmtab, spmij, SCspm, dijmn, dijmm]
Clear[kx, ky]
Hsc = H*Exp[(*I*Re[kx]*Nx*Phasex*) I*Re[ky]*Ny*Phasey]; 
(*hudsc = hud*Exp[I*Re[kx]*Nx*Phasexud + I*Re[ky]*Ny*Phaseyud]; 
hdusc = hdu*Exp[I*Re[kx]*Nx*Phasexdu + I*Re[ky]*Ny*Phaseydu]; 
hdnsc = hdn*Exp[I*Re[kx]*Nx*Phasexdn + I*Re[ky]*Ny*Phaseydn]; 
H=ArrayFlatten[{{hupsc,hudsc},{hdusc,hdnsc}}];
H == ConjugateTranspose[H]; *)
Clear[hup,hud,hdu,hdn,rhup,ihup,rhud,ihud,rhdu,ihdu,rhdn,ihdn]; 



(*kkmx = Ns*Nx; *)
kkmy = Ns*Ny; 


(*kT=0.0035;
ft[x_]:=(*1./(Exp[Re[x]/kT]+1.)-1./(Exp[Re[x]/kT]+1.)^2*)If[Re[x]>3,0.0,
If[Re[x]<-3,0.0,
1./(Exp[Re[x]/kT]+1.)-1./(Exp[Re[x]/kT]+1.)^2]];*)

Print["starting ldos calc"]

 
 CloseKernels[];
 LaunchKernels[kernum];
ldos=Total[(*Flatten[*)ParallelTable[
  ES = Eigensystem[Hsc];
   elist = Table[ES[[1,l]], {l, 1, Nband*NN}];
    ulist = Table[Abs[ES[[2,l,ival]]]^2, {l, 1, Nband*NN}]; 
      (*vlist = Table[Abs[ES[[2,l,i + Nband/2*NN*2]]]^2, {l, 1, Nband/2*4*NN}, {i, 388+1, 388+(nsite - 1)*Nband/2 + Nband/2}];*) 
      
      ldossc = Table[Im[Total[Table[ulist[[l,1 ;; All]]*(1/(-wmax + wmax*2*(w/nw) - elist[[l]] + I*\[Eta]))(* + 
      vlist[[l,1 ;; All]]*(1/(-wmax + 2*wmax*(w/nw) + elist[[l]] + I*\[Eta]))*), {l, 1, Nband*NN}]]], {w, 0, nw}]; 
      
      ldossc(*,{kx,0,(Ns - 1.)*2*(Pi/kkmx),2*(Pi/kkmx)}*),{ky,0,(Ns - 1.)*2*(Pi/kkmy),2*(Pi/kkmy)}](*,1]*)];
      Print["done"];
      
      (*ldossc=Table[Total[Table[ulist[[l,1 ;; All]]*ft[-wmax + wmax*2*(w/nw) - elist[[l]]] (*+
vlist[[l,1 ;; All]]*ft[-wmax + 2*wmax*(w/nw) + elist[[l]]]*), {l, 1,Nband*NN}]], {w, 0, nw}];
ldossc,{kx,0,(Ns - 1.)*2*(Pi/kkmx),2*(Pi/kkmx)},{ky,0,(Ns - 1.)*2*(Pi/kkmy),2*(Pi/kkmy)}],1]];

Print["done"];*)
      
    
  (*Export["ldos_wte2_cdw.dat",(* ldos*)Table[{-wmax + wmax*2*((n-1)/nw),(-Pi^(-1))*ldos[[n,i]]/Ns^2}, {n, 1, nw+1}, {i, 1, Length[ival]}]]
Print["done"];*)
fl = StringJoin["yedge_mu",ToString[mu],"_cdw",ToString[cdw],"_wmax",ToString[wmax],"Ny",ToString[Ny],"Nx",ToString[Nx],
"eta",ToString[\[Eta]],"Ns",ToString[Ns]];
DeleteFile[StringJoin[fl,"Ldos_wte2_cdw.dat"]];
str1 = OpenAppend[StringJoin[fl,"Ldos_wte2_cdw.dat"], FormatType -> StandardForm];
For[n=1,n<=nw+1,n++,
For[i=1,i<=Length[ival],i++,
Write[str1,-wmax + wmax*2*((n-1)/nw),"  ", (-1/(Pi))*ldos[[n,i]]/Ns];
]
];
Close[str1];
Exit[];
     






