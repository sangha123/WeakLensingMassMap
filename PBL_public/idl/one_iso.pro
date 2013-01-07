pro one_iso,nx=nx,theta_e=theta_e,theta_c=theta_c

; Run this first.  This generates the cluster
; For our fiducial value, we've been taking:
; theta_e=0.2, theta_c=0.1 for both Softened isothermals

if not keyword_set(nx) then nx=100
if not keyword_set(theta_e) then theta_e=0.25
if not keyword_set(theta_c) then theta_c=0.12

theta_e=theta_e*nx
theta_c=theta_c*nx

xc=nx/2.
yc=nx/2.

kappa=dblarr(nx,nx)*0
psi=dblarr(nx,nx)*0
alpha1=dblarr(nx,nx)*0
alpha2=dblarr(nx,nx)*0
gamma1=dblarr(nx,nx)*0
gamma2=dblarr(nx,nx)*0
beta1=dblarr(nx,nx)*0
beta2=dblarr(nx,nx)*0
theta_softsq=dblarr(nx,nx)*0
thetasq=dblarr(nx,nx)*0

for i=0,nx-1 do begin
   for j=0,nx-1 do begin
      thetasq[i,j]=((i-xc)^2+(j-yc)^2) 
      theta_softsq[i,j]=thetasq[i,j]+theta_c^2
      kappa[i,j]=0.5*theta_e*(thetasq[i,j]+2*theta_c^2)/theta_softsq[i,j]^1.5
      psi[i,j]=theta_e*sqrt(theta_softsq[i,j])
      alpha1[i,j]=(theta_e/(sqrt(theta_softsq[i,j])))*(i-xc)
      alpha2[i,j]=(theta_e/(sqrt(theta_softsq[i,j])))*(j-yc)
      gamma1[i,j]=(theta_e*((j-yc)^2-(i-xc)^2))/(2.*theta_softsq[i,j]^1.5)
      gamma2[i,j]=-(theta_e*(j-yc)*(i-xc)/(theta_softsq[i,j]^1.5))
      beta1[i,j]=((i-nx/2.)-alpha1[i,j])
      beta2[i,j]=((j-nx/2.)-alpha2[i,j])
   endfor
endfor
flux=1./((1.-kappa)^2-(gamma1^2+gamma2^2)) 
cluster={ theta_e:theta_e,  $
          theta_c:theta_c,  $
          xc:xc,            $
          yc:yc,            $
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

save,file="cluster_1peak.sav",cluster
end
 
