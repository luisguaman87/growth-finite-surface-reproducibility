# Growth on a Finite Surface — reproducible numerical illustrations
# Guamán-Lazo, Luis | Universidad Técnica de Manabí
# Base R only. From repository root: source('code/Reproducible_model_code.R')
# Outputs are written in ./finite_surface_outputs/
# NOTE: equilibrium theorems concern a REDUCED constant-saving system,
# not existence/uniqueness of the infinite-horizon planner optimum.

set.seed(20260928)
# Run from the repository root. No machine-specific working directory is required.
outdir <- 'finite_surface_outputs'
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
p <- list(Tbar=1, fbar=.18, AF=.65, xi=.50, hbar=.12, AH=.80,
          A=.90, alpha=.35, beta=.25, gamma=.20, s=.22, delta=.05,
          rhoN=.08, Nbar=1, chi=.015, zeta=.005, Ncrit=.20,
          LN_min=.12)
# Explicit parameter alias: prevents undefined defaults in helper functions.
PARAMS <- p
stopifnot(p$Tbar>0,p$fbar>0,p$AF>0,p$xi>0,p$hbar>0,p$AH>0,
          p$A>0,p$alpha>0,p$alpha<1,p$beta>0,p$gamma>=0,
          p$s>0,p$s<1,p$delta>0,p$rhoN>0,p$Nbar>0,
          p$chi>0,p$zeta>0,p$Ncrit>0,p$Ncrit<p$Nbar,
          p$LN_min>=0,p$LN_min<p$Tbar)

land <- function(P,N,LN,p=PARAMS) {
  if(any(!is.finite(c(P,N,LN))) || P<=0 || N<=0 || LN<p$LN_min || LN>p$Tbar)
    stop('Invalid population, ecological stock or conservation allocation')
  LA <- p$fbar*P/(p$AF*N^p$xi)
  LH <- p$hbar*P/p$AH
  LY <- p$Tbar-LA-LH-LN
  c(LA=LA,LH=LH,LY=LY,LN=LN)
}
D <- function(P,N,LN,p=PARAMS) unname(land(P,N,LN,p)['LY'])
Y <- function(K,P,N,LN,p=PARAMS) {
  if(!is.finite(K) || K<=0 || D(P,N,LN,p)<=0) return(NA_real_)
  p$A*K^p$alpha*D(P,N,LN,p)^p$beta*N^p$gamma
}
dD_dP <- function(P,N,p=PARAMS) -p$fbar/(p$AF*N^p$xi)-p$hbar/p$AH
dD_dN <- function(P,N,p=PARAMS) p$xi*p$fbar*P/(p$AF*N^(p$xi+1))
d2D_dN2 <- function(P,N,p=PARAMS) -p$xi*(p$xi+1)*p$fbar*P/(p$AF*N^(p$xi+2))
d2D_dPdN <- function(P,N,p=PARAMS) p$xi*p$fbar/(p$AF*N^(p$xi+1))
dY_dK <- function(K,P,N,LN,p=PARAMS) p$alpha*Y(K,P,N,LN,p)/K
# IMPORTANT: planner partial holds LY fixed; reduced derivative includes D'(N).
dY_dN_partial <- function(K,P,N,LN,p=PARAMS) p$gamma*Y(K,P,N,LN,p)/N
dY_dN_total <- function(K,P,N,LN,p=PARAMS) {
  y <- Y(K,P,N,LN,p)
  y*(p$gamma/N+p$beta*dD_dN(P,N,p)/D(P,N,LN,p))
}
Pmax <- function(N,LN=p$LN_min,p=PARAMS) (p$Tbar-LN)/(p$fbar/(p$AF*N^p$xi)+p$hbar/p$AH)
N_feasibility <- function(P,LN,p=PARAMS) {
  B <- p$Tbar-p$hbar*P/p$AH-LN
  if(B<=0) return(NA_real_)
  (p$fbar*P/(p$AF*B))^(1/p$xi)
}
cdiff <- function(fun,x,h=1e-5) (fun(x+h)-fun(x-h))/(2*h)
P0 <- .75; N0 <- .80; K0 <- 1.50; LN0 <- .22
stopifnot(LN0>=p$LN_min,D(P0,N0,LN0,p)>0)
checks <- data.frame(
  derivative=c('D_P','D_N','D_NN','D_PN','D_LN','Y_K','Y_N_total','Y_N_partial'),
  analytical=c(dD_dP(P0,N0,p),dD_dN(P0,N0,p),d2D_dN2(P0,N0,p),
               d2D_dPdN(P0,N0,p),-1,dY_dK(K0,P0,N0,LN0,p),
               dY_dN_total(K0,P0,N0,LN0,p),dY_dN_partial(K0,P0,N0,LN0,p)),
  numerical=c(cdiff(function(x) D(x,N0,LN0,p),P0),
              cdiff(function(x) D(P0,x,LN0,p),N0),
              cdiff(function(x) dD_dN(P0,x,p),N0),
              cdiff(function(x) dD_dP(P0,x,p),N0),
              cdiff(function(x) D(P0,N0,x,p),LN0),
              cdiff(function(x) Y(x,P0,N0,LN0,p),K0),
              cdiff(function(x) Y(K0,P0,x,LN0,p),N0),
              cdiff(function(x) p$A*K0^p$alpha*D(P0,N0,LN0,p)^p$beta*x^p$gamma,N0)))
checks$absolute_error <- abs(checks$analytical-checks$numerical)
write.csv(checks,file.path(outdir,'derivative_checks.csv'),row.names=FALSE)
stopifnot(all(is.finite(checks$absolute_error)),max(checks$absolute_error)<1e-5)

# Proposition 3: analytic Nf, strict boundary condition and unique root.
solve_reduced <- function(P,LN,p=PARAMS) {
  nf <- N_feasibility(P,LN,p)
  if(!is.finite(nf) || nf>=p$Nbar)
    return(list(status='No nonempty feasible ecological interval',Nf=nf))
  margin <- p$rhoN*LN*(1-nf/p$Nbar)-p$zeta*p$fbar*P
  if(margin<=0)
    return(list(status='Existence theorem boundary condition fails',Nf=nf,margin=margin))
  Knull <- function(N) ((p$s*p$A/p$delta)*D(P,N,LN,p)^p$beta*N^p$gamma)^(1/(1-p$alpha))
  G <- function(N) p$rhoN*LN*(1-N/p$Nbar)-p$chi*Y(Knull(N),P,N,LN,p)-p$zeta*p$fbar*P
  # G is evaluated at Nf by its analytic right-hand limit, avoiding D=0.
  ns <- uniroot(G,c(nf+1e-10*(p$Nbar-nf),p$Nbar),tol=1e-12)$root
  ks <- Knull(ns); ys <- Y(ks,P,ns,LN,p)
  YK <- dY_dK(ks,P,ns,LN,p); YN <- dY_dN_total(ks,P,ns,LN,p)
  J <- matrix(c(p$s*YK-p$delta,-p$chi*YK,p$s*YN,
                -p$rhoN*LN/p$Nbar-p$chi*YN),2,2,
              dimnames=list(c('Kdot','Ndot'),c('K','N')))
  tr <- sum(diag(J)); determinant <- det(J); eigenvalues <- eigen(J)$values
  stopifnot(tr<0,determinant>0,all(Re(eigenvalues)<0),
            abs(p$s*ys-p$delta*ks)<1e-7,abs(G(ns))<1e-7)
  list(status=if(ns>p$Ncrit) 'Interior and ecologically admissible' else
         'Interior mathematically; violates Ncrit',Nf=nf,margin=margin,
       Nstar=ns,Kstar=ks,Ystar=ys,J=J,trace=tr,determinant=determinant,
       eigenvalues=eigenvalues,Knull=Knull,G=G)
}
sol <- solve_reduced(P0,LN0,p)
if(is.null(sol$Nstar)) stop(sol$status)
write.csv(data.frame(Nf=sol$Nf,Nstar=sol$Nstar,Kstar=sol$Kstar,
                     Ystar=sol$Ystar,trace=sol$trace,determinant=sol$determinant,
                     admissible=sol$Nstar>p$Ncrit),
          file.path(outdir,'equilibrium_v6.csv'),row.names=FALSE)
write.csv(sol$J,file.path(outdir,'jacobian_v6.csv'))

# RK4 with no artificial floor or clipping: stop at first detected infeasible step (grid approximation).
rhs <- function(z,LN,P,p=PARAMS) {
  K <- z[1]; N <- z[2]
  if(!is.finite(K) || !is.finite(N) || K<=0 || N<=0 || N>p$Nbar ||
     D(P,N,LN,p)<=0) return(c(NA_real_,NA_real_))
  yy <- Y(K,P,N,LN,p)
  c(p$s*yy-p$delta*K,
    p$rhoN*LN*(1-N/p$Nbar)-p$chi*yy-p$zeta*p$fbar*P)
}
simulate_path <- function(LN,T=150,dt=.05,P=P0,p=PARAMS,Kinit=1,Ninit=.85) {
  tt <- seq(0,T,by=dt); zz <- matrix(NA_real_,length(tt),2)
  zz[1,] <- c(Kinit,Ninit); status <- 'completed'; last <- length(tt)
  for(i in 2:length(tt)) {
    z <- zz[i-1,]; a <- rhs(z,LN,P,p)
    b <- rhs(z+dt*a/2,LN,P,p)
    c <- rhs(z+dt*b/2,LN,P,p)
    d <- rhs(z+dt*c,LN,P,p)
    zn <- z+dt*(a+2*b+2*c+d)/6
    if(any(!is.finite(c(a,b,c,d,zn))) || zn[1]<=0 || zn[2]<=0 ||
       zn[2]>p$Nbar || D(P,zn[2],LN,p)<=0 || zn[2]<p$Ncrit) {
      status <- 'stopped: infeasible or ecological floor reached'
      last <- i-1; break
    }
    zz[i,] <- zn
  }
  zz <- zz[seq_len(last),,drop=FALSE]
  zN <- zz[,2]
  yy <- vapply(seq_len(last),function(i) Y(zz[i,1],P,zN[i],LN,p),numeric(1))
  LL <- t(vapply(zN,function(n) land(P,n,LN,p),numeric(4)))
  df <- data.frame(time=tt[seq_len(last)],K=zz[,1],N=zN,Y=yy,
                   LA=LL[,1],LH=LL[,2],LY=LL[,3],LN=LL[,4],
                   scenario_LN=LN)
  stopifnot(max(abs(rowSums(LL)-p$Tbar))<1e-10,all(df$LY>0))
  attr(df,'status') <- status
  df
}
scenarios <- c(.15,.22,.30)
paths <- lapply(scenarios, function(LN) simulate_path(LN=LN, P=P0, p=PARAMS))
allpaths <- do.call(rbind,paths)
write.csv(allpaths,file.path(outdir,'transition_paths.csv'),row.names=FALSE)
write.csv(data.frame(LN=scenarios,status=vapply(paths,function(x) attr(x,'status'),''),
                     last_time=vapply(paths,function(x) tail(x$time,1),0.0)),
          file.path(outdir,'transition_status.csv'),row.names=FALSE)

# Figures: reproducible 600 dpi TIFF + PNG + vector PDF from the SAME plotting function.
colors <- c('#284B63','#4C956C','#E4A44B','#B56B45')
export_plot <- function(stem,draw,width=7.5,height=5.2) {
  for(fmt in c('png','tiff','pdf')) {
    fn <- file.path(outdir,paste0(stem,'.',fmt))
    if(fmt=='png') png(fn,width=width,height=height,units='in',res=600)
    if(fmt=='tiff') tiff(fn,width=width,height=height,units='in',res=600,compression='lzw')
    if(fmt=='pdf') pdf(fn,width=width,height=height,onefile=FALSE)
    draw(); dev.off()
  }
}
# Publication graphics: compact axis/legend sizes, consistent with manuscript order.
figure_phase <- function() {
  nv <- seq(sol$Nf+.002,p$Nbar-.002,length.out=600)
  kn <- vapply(nv,sol$Knull,numeric(1))
  nn <- vapply(nv,function(n) {
    regen <- p$rhoN*LN0*(1-n/p$Nbar)-p$zeta*p$fbar*P0
    if(regen<=0) return(NA_real_)
    (regen/(p$chi*p$A*D(P0,n,LN0,p)^p$beta*n^p$gamma))^(1/p$alpha)
  },numeric(1))
  par(mar=c(4.6,4.9,1.2,1.2),las=1,bty='l')
  plot(nv,kn,type='l',lwd=2.1,col=colors[1],ylim=c(0,max(8,1.4*sol$Kstar)),
       xlab=expression('Functional natural capital, '*italic(N)*' (normalized)'),
       ylab=expression('Produced capital, '*italic(K)*' (normalized)'),cex.lab=1.00,cex.axis=.90)
  lines(nv,nn,lwd=2.1,lty=2,col=colors[2]); points(sol$Nstar,sol$Kstar,pch=19,col=colors[3],cex=1.2)
  text(sol$Nstar,sol$Kstar,labels=expression(group('(',list(N^'*',K^'*'),')')),
       pos=4,offset=.8,cex=.94)
  abline(v=p$Ncrit,lty=3,lwd=1.2,col='grey50')
  legend('topright',bty='n',lwd=c(2.1,2.1,NA,1.2),lty=c(1,2,NA,3),
         pch=c(NA,NA,19,NA),col=c(colors[1:3],'grey50'),
         legend=as.expression(list(bquote(dot(K)==0),bquote(dot(N)==0),
                                   'Steady-state equilibrium',bquote(N[crit]))),cex=.88)
}
figure_land <- function() {
  Pgrid <- seq(.45,1.30,length.out=250)
  MM <- t(vapply(Pgrid,function(P) land(P,N0,LN0,p),numeric(4)))
  stopifnot(all(MM[,'LY']>0)); MM <- MM/p$Tbar
  cum <- t(apply(MM,1,cumsum))
  par(mar=c(7.3,4.7,1.2,1.2),las=1,bty='l')
  plot(Pgrid,Pgrid*0,type='n',ylim=c(0,1),xlab='',
       ylab='Share of total land',cex.lab=1.00,cex.axis=.90)
  for(j in seq_len(4)) {
    low <- if(j==1) rep(0,length(Pgrid)) else cum[,j-1]
    polygon(c(Pgrid,rev(Pgrid)),c(low,rev(cum[,j])),col=colors[j],border=NA)
  }
  mtext(expression('Population, '*italic(P)*' (normalized)'),side=1,line=2.6,cex=1.00)
  legend('bottom',inset=c(0,-.45),xpd=NA,ncol=2,bty='n',
         legend=expression('Agriculture ('*L[A]*')','Housing ('*L[H]*')',
                           'Market production ('*L[Y]*')','Natural land ('*L[N]*')'),
         fill=colors,cex=.88,x.intersp=.8,y.intersp=1.12)
}
figure_paths <- function(field='Y') {
  par(mar=c(4.5,4.8,1.2,1.2),las=1,bty='l')
  yr <- range(unlist(lapply(paths,function(z) z[[field]])),finite=TRUE)
  yl <- switch(field,Y=expression('Market output, '*italic(Y)*' (normalized)'),
               N=expression('Functional natural capital, '*italic(N)*' (normalized)'),
               LY=expression('Market land, '*L[Y]*' (normalized)'),field)
  plot(paths[[1]]$time,paths[[1]][[field]],type='l',lwd=2,col=colors[1],
       ylim=yr,xlab='Time (model units)',ylab=yl,cex.lab=1.00,cex.axis=.90)
  for(i in 2:length(paths)) lines(paths[[i]]$time,paths[[i]][[field]],
                                  lwd=2,col=colors[i],lty=i)
  legend('topright',bty='n',lwd=2,lty=seq_along(paths),col=colors[seq_along(paths)],
         legend=as.expression(lapply(scenarios,function(z) bquote(L[N]==.(z)))),cex=.88)
}
figure_release <- function() {
  nv <- seq(max(p$Ncrit,.22),p$Nbar,length.out=250)
  ag <- vapply(nv,function(n) unname(land(P0,n,LN0,p)['LA']),numeric(1))/p$Tbar
  pr <- vapply(nv,function(n) D(P0,n,LN0,p),numeric(1))/p$Tbar
  par(mar=c(4.5,4.8,1.2,1.2),las=1,bty='l')
  plot(nv,ag,type='l',lwd=2,col=colors[1],ylim=range(c(ag,pr)),
       xlab=expression('Functional natural capital, '*italic(N)*' (normalized)'),
       ylab='Share of total land',cex.lab=1.00,cex.axis=.90)
  lines(nv,pr,lwd=2,col=colors[2],lty=2)
  legend('right',bty='n',lwd=2,lty=c(1,2),col=colors[1:2],
         legend=expression('Minimum agricultural land ('*L[A]^min*')',
                           'Residual market land ('*L[Y]*')'),
         cex=.88,x.intersp=.8,y.intersp=1.12)
}
# Manuscript numbering: phase, population displacement, output paths, ecological release.
export_plot('Figure_1_phase_diagram',figure_phase)
export_plot('Figure_2_land_displacement',figure_land)
export_plot('Figure_3_output_paths',function() figure_paths('Y'))
export_plot('Figure_4_ecological_land_release',figure_release)
# Supplementary diagnostics, not numbered as main-text figures.
export_plot('Supplement_Natural_capital_paths',function() figure_paths('N'))
export_plot('Supplement_Productive_land_paths',function() figure_paths('LY'))

# Sensitivity: distinguish theoretical existence, ecological admissibility and failure.
grid <- expand.grid(LN=c(.15,.20,.25,.30),rhoN=c(.05,.08,.11),
                    chi=c(.010,.015,.020),xi=c(.25,.50,.75))
sensitivity <- do.call(rbind,lapply(seq_len(nrow(grid)),function(i) {
  pp <- p; pp$rhoN <- grid$rhoN[i]; pp$chi <- grid$chi[i]; pp$xi <- grid$xi[i]
  ss <- solve_reduced(P0,grid$LN[i],pp)
  data.frame(grid[i,],status=ss$status,Nf=ss$Nf,
             Nstar=if(is.null(ss$Nstar)) NA_real_ else ss$Nstar,
             Kstar=if(is.null(ss$Kstar)) NA_real_ else ss$Kstar)
}))
write.csv(sensitivity,file.path(outdir,'sensitivity.csv'),row.names=FALSE)
cat('Baseline mathematical checks passed. Output directory:',normalizePath(outdir),'\n')
cat(sprintf('Nf=%.6f; N*=%.6f; K*=%.6f; Y*=%.6f; ecological admissible=%s\n',
            sol$Nf,sol$Nstar,sol$Kstar,sol$Ystar,sol$Nstar>p$Ncrit))

# MODEL A vs MODEL B: theoretical robustness, gamma=0 (not a planner optimum).
model_comparison <- do.call(rbind,lapply(c('A_direct_and_indirect','B_indirect_only'),function(model) {
  pp <- p; if(model=='B_indirect_only') pp$gamma <- 0
  ss <- solve_reduced(P0,LN0,pp)
  if(is.null(ss$Nstar)) return(data.frame(model=model,gamma=pp$gamma,status=ss$status,Nstar=NA_real_,Kstar=NA_real_,Ystar=NA_real_,trace=NA_real_,determinant=NA_real_))
  data.frame(model=model,gamma=pp$gamma,status=ss$status,Nstar=ss$Nstar,Kstar=ss$Kstar,Ystar=ss$Ystar,trace=ss$trace,determinant=ss$determinant)
}))
write.csv(model_comparison,file.path(outdir,'model_A_vs_B_gamma0.csv'),row.names=FALSE)
print(model_comparison)
# Dulac sign at benchmark interior grid: analytic sign check, not global-stability proof.
dulac_check <- function(K,N,P,LN,p=PARAMS) {
  yy <- Y(K,P,N,LN,p); yn <- dY_dN_total(K,P,N,LN,p)
  p$s*(p$alpha-1)*yy/K^2-(p$rhoN*LN/p$Nbar+p$chi*yn)/K
}
dulac_grid <- expand.grid(K=seq(.1,8,length.out=30),N=seq(sol$Nf+.01,p$Nbar-.01,length.out=30))
dulac_grid$divergence <- mapply(function(K,N) dulac_check(K,N,P0,LN0,p),dulac_grid$K,dulac_grid$N)
stopifnot(all(is.finite(dulac_grid$divergence)),all(dulac_grid$divergence<0))
write.csv(dulac_grid,file.path(outdir,'dulac_negative_divergence_grid.csv'),row.names=FALSE)

# --- Exogenous demographic sensitivity (not an autonomous stability theorem) ---
# The benchmark population P0 remains fixed in all existing equilibrium results.
# Population varies ONLY in this additional, explicitly nonautonomous exercise.
simulate_demography <- function(gP,LN=LN0,T=100,dt=.05,Pinit=P0,Kinit=1,Ninit=.85,p=PARAMS) {
  stopifnot(gP>=0,LN>=p$LN_min,dt>0,Pinit>0)
  times <- seq(0,T,by=dt); z <- c(Kinit,Ninit)
  data <- data.frame(time=numeric(0),P=numeric(0),K=numeric(0),N=numeric(0),Y=numeric(0),LY=numeric(0),gP=numeric(0))
  status <- 'completed'
  pmax_cap <- Pmax(p$Nbar,LN,p)
  for(i in seq_along(times)) {
    t <- times[i]; pop <- Pinit*exp(gP*t)
    if(!all(is.finite(z)) || z[1]<=0 || z[2]<p$Ncrit || z[2]>p$Nbar || D(pop,z[2],LN,p)<=0) {
      status <- 'stopped: feasibility boundary'; break
    }
    data <- rbind(data,data.frame(time=t,P=pop,K=z[1],N=z[2],Y=Y(z[1],pop,z[2],LN,p),LY=D(pop,z[2],LN,p),gP=gP))
    if(i==length(times)) break
    f <- function(tt,zz) rhs(zz,LN,Pinit*exp(gP*tt),p)
    a <- f(t,z); b <- f(t+dt/2,z+dt*a/2)
    c <- f(t+dt/2,z+dt*b/2); d <- f(t+dt,z+dt*c)
    zn <- z+dt*(a+2*b+2*c+d)/6
    if(any(!is.finite(c(a,b,c,d,zn))) || zn[1]<=0 || zn[2]<p$Ncrit || zn[2]>p$Nbar || D(Pinit*exp(gP*(t+dt)),zn[2],LN,p)<=0) {
      status <- 'stopped: feasibility boundary'; break
    }
    z <- zn
  }
  # Last feasible grid time is not an exact root-located event time.
  data$pressure <- data$P/pmax_cap
  attr(data,'status') <- status
  data
}
demo_rates <- c(0,.005,.01)
demo_paths <- lapply(demo_rates,function(g) simulate_demography(gP=g,LN=LN0,p=PARAMS))
write.csv(do.call(rbind,demo_paths),file.path(outdir,'demographic_sensitivity.csv'),row.names=FALSE)
write.csv(data.frame(gP=demo_rates,status=vapply(demo_paths,function(z) attr(z,'status'),''),
                     last_feasible_time=vapply(demo_paths,function(z) tail(z$time,1),0.0),
                     theoretical_upper_bound=vapply(demo_rates,function(g) if(g==0) Inf else log(Pmax(p$Nbar,LN0,p)/P0)/g,0.0)),
          file.path(outdir,'demographic_sensitivity_status.csv'),row.names=FALSE)
figure_demographic <- function() {
  old <- par(no.readonly=TRUE); on.exit(par(old))
  par(mfrow=c(1,2),mar=c(4.4,4.6,1.5,1),las=1,bty='l')
  yr <- range(unlist(lapply(demo_paths,function(z) z$pressure)),finite=TRUE)
  plot(NA,xlim=c(0,100),ylim=range(c(yr,1)),xlab='Time (model units)',
       ylab=expression('Population pressure, '*P/P[max]),cex.lab=1.00,cex.axis=.90)
  abline(h=1,col='grey60',lty=3)
  for(i in seq_along(demo_paths)) lines(demo_paths[[i]]$time,demo_paths[[i]]$pressure,col=colors[i],lty=i,lwd=2)
  legend('topleft',legend=expression(g[P]==0,g[P]==0.005,g[P]==0.01),col=colors[1:3],lty=1:3,lwd=2,bty='n',cex=.85)
  yr <- range(unlist(lapply(demo_paths,function(z) z$Y)),finite=TRUE)
  plot(NA,xlim=c(0,100),ylim=yr,xlab='Time (model units)',
       ylab=expression('Non-agricultural market output, '*Y*' (normalized)'),cex.lab=1.00,cex.axis=.90)
  for(i in seq_along(demo_paths)) lines(demo_paths[[i]]$time,demo_paths[[i]]$Y,col=colors[i],lty=i,lwd=2)
  legend('bottomleft',legend=expression(g[P]==0,g[P]==0.005,g[P]==0.01),col=colors[1:3],lty=1:3,lwd=2,bty='n',cex=.85)
}
export_plot('Figure_5_demographic_sensitivity',figure_demographic,width=10.8,height=4.4)
cat('Demographic sensitivity exported; benchmark stability proof applies only to gP=0.\n')
