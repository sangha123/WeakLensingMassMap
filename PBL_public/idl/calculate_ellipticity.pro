pro calculate_ellipticity,x,y,es1,es2,e1,e2,infile=infile
;;;es1,es2 are the two components of intrinsic ellipticity
;;;infile= .sav file having the profile
;;; This code calculates the induced elliticity at a given x,y

if not keyword_set(infile) then infile='cluster.sav'

restore,infile
  
g1=cluster.gamma1[x,y]/(1-cluster.kappa[x,y])
g2=cluster.gamma2[x,y]/(1-cluster.kappa[x,y])

g=sqrt(g1^2+g2^2)

if g lt 1 then begin
   A=es1+g1
   B=es2+g2
   C=1+g1*es1+g2*es2
   D=g1*es2-g2*es1              
endif else begin
   A=1+g1*es1+g2*es2            
   B=g2*es1-g1*es2
   C=es1+g1
   D=-es2-g2
endelse
    e1=(A*C+B*D)/(C*C+D*D)
    e2=(B*C-A*D)/(C*C+D*D)
end
