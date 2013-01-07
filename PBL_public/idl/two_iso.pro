pro two_iso,nx=nx,theta_e1=theta_e1,theta_c1=theta_c1,theta_e2=theta_e2,$
                 theta_c2=theta_c2,xc1=xc1,yc1=yc1,xc2=xc2,yc2=yc2

; Run this first.  This generates the cluster
; For our fiducial value, we've been taking:
; theta_e=0.2, theta_c=0.1 for both Softened isothermals


if not keyword_set(nx) then nx=100
if not keyword_set(theta_e1) then theta_e1=0.2
if not keyword_set(theta_e2) then theta_e2=0.2
if not keyword_set(theta_c1) then theta_c1=0.1
if not keyword_set(theta_c2) then theta_c2=0.1
if not keyword_set(xc1) then xc1=0.35
if not keyword_set(yc1) then yc1=0.7
if not keyword_set(xc2) then xc2=0.7
if not keyword_set(yc2) then yc2=0.35


theta_e1=theta_e1*nx
theta_e2=theta_e2*nx
theta_c1=theta_c1*nx
theta_c2=theta_c2*nx
xc1=xc1*nx
yc1=yc1*nx
xc2=xc2*nx
yc2=yc2*nx


kappa1=dblarr(nx,nx)*0
psi1=dblarr(nx,nx)*0
alpha1x=dblarr(nx,nx)*0
alpha1y=dblarr(nx,nx)*0
gamma1x=dblarr(nx,nx)*0
gamma1y=dblarr(nx,nx)*0
theta_softsq1=dblarr(nx,nx)*0
thetasq1=dblarr(nx,nx)*0
kappa2=dblarr(nx,nx)*0
psi2=dblarr(nx,nx)*0
alpha2x=dblarr(nx,nx)*0
alpha2y=dblarr(nx,nx)*0
gamma2x=dblarr(nx,nx)*0
gamma2y=dblarr(nx,nx)*0

for i=0,nx-1 do begin
   for j=0,nx-1 do begin
      thetasq1[i,j]=((i-xc1)^2+(j-yc1)^2) 
      theta_softsq1[i,j]=thetasq1[i,j]+theta_c1^2
      kappa1[i,j]=0.5*theta_e1*(thetasq1[i,j]+2*theta_c1^2)/theta_softsq1[i,j]^1.5
      psi1[i,j]=theta_e1*sqrt(theta_softsq1[i,j])
      alpha1x[i,j]=(theta_e1/(sqrt(theta_softsq1[i,j])))*(i-xc1)
      alpha1y[i,j]=(theta_e1/(sqrt(theta_softsq1[i,j])))*(j-yc1)
      gamma1x[i,j]=(theta_e1*((j-yc1)^2-(i-xc1)^2))/(2.*theta_softsq1[i,j]^1.5)
      gamma1y[i,j]=-(theta_e1*(j-yc1)*(i-xc1)/(theta_softsq1[i,j]^1.5))

      thetasq2=((i-xc2)^2+(j-yc2)^2)
      theta_softsq2=thetasq2+theta_c2^2
      kappa2[i,j]=0.5*theta_e2*(thetasq2+2*theta_c2^2)/theta_softsq2^1.5
      psi2[i,j]=theta_e2*sqrt(theta_softsq2)
      alpha2x[i,j]=(theta_e2/(sqrt(theta_softsq2)))*(i-xc2)
      alpha2y[i,j]=(theta_e2/(sqrt(theta_softsq2)))*(j-yc2)
      gamma2x[i,j]=(theta_e2*((j-yc2)^2-(i-xc2)^2))/(2.*theta_softsq2^1.5)
      gamma2y[i,j]=-(theta_e2*(j-yc2)*(i-xc2)/(theta_softsq2^1.5))
   endfor
endfor



psi=psi1+psi2
kappa=kappa1+kappa2
alpha1=alpha1x+alpha2x
alpha2=alpha1y+alpha2y
gamma1=gamma1x+gamma2x
gamma2=gamma1y+gamma2y

beta1=dblarr(nx,nx)
beta2=dblarr(nx,nx)

for i=0,nx-1 do begin
   for j=0,nx-1 do begin
      beta1[i,j]=((i-nx/2.)-alpha1[i,j])
      beta2[i,j]=((j-nx/2.)-alpha2[i,j])
   endfor
endfor

flux=1./((1.-kappa)^2-(gamma1^2+gamma2^2)) 

cluster={ theta_e1:theta_e1,  $
          theta_c1:theta_c1,  $
          theta_e2:theta_e2,  $
          theta_c2:theta_c2,  $
          xc1:xc1,            $
          yc1:yc1,            $
          xc2:xc2,            $
          yc2:yc2,            $
          kappa:kappa,        $
          psi:psi,            $
          alpha1:alpha1,      $ ;;;the deflection field        
          alpha2:alpha2,      $
          beta1:beta1,        $ ;;;the image field
          beta2:beta2,        $
          gamma1:gamma1,      $
          gamma2:gamma2,      $
          flux:flux,          $
          nx:nx}

save,file="cluster.sav",cluster

end
