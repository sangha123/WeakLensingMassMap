pro find_images,ngal=ngal,outfile=outfile,sig_e=sig_e,infile=infile

;;;ngal:number of background galaxies.

if not keyword_set(outfile) then outfile="singmult.dat"
if not keyword_set(ngal) then ngal=1000
if not keyword_set(sig_e) then sig_e=0.1
if not keyword_set(infile) then infile='cluster.sav'
;;;WARNING: THis needs to be taken care of... if two peak change name of .sav
;;; file and function name. It also needS to be changed in calculate_ellipticity...

print,systime()

restore,infile,/verbose

nx=cluster.nx

xc=(nx*randomu(seed,ngal))-nx/2.
yc=(nx*randomu(seed,ngal))-nx/2.

es1=sig_e*randomn(5,ngal)
es2=sig_e*randomn(7,ngal)

delta=1.
;; ;;reduced shear
;; gamma_r=complex(cluster.gamma1,cluster.gamma2)/(1-cluster.kappa)
close,1
openw,1,outfile
num=0
nm=0

g1=cluster.gamma1/(1-cluster.kappa)
g2=cluster.gamma2/(1-cluster.kappa)

for i=0,ngal-1 do begin
   
   ;read,idum
    idx=where(cluster.beta1 gt xc[i]-delta and cluster.beta1 lt xc[i]+delta and cluster.beta2 gt yc[i]-delta and cluster.beta2 lt yc[i]+delta,nloc)
    ;print,xc[i],yc[i],nloc
    if nloc gt 0 then begin
       nm=nm+1
       i_pos=(idx mod nx)
       j_pos=(idx/nx)
       ;plot,i_pos,j_pos,psym=4,xrange=[0,nx],yrange=[0,nx]
       
       find_clump,nx,i_pos,j_pos,xc_img,yc_img,n
       
       
       mult=0.
       id=where((xc_img-nx/2)^2+(yc_img-nx/2)^2 gt 25,n)
       if n gt 0 then begin 
          xc_img=xc_img[id]
          yc_img=yc_img[id]
       endif
       if n gt 1 then begin
          size_img=(size(xc_img))[1]
          
          del_x=(xc_img-cluster.alpha1[nint(xc_img),nint(yc_img)]-nx/2.)
          del_y=(yc_img-cluster.alpha2[nint(xc_img),nint(yc_img)]-nx/2.)
          del2_x=abs(del_x-xc[i])
          del2_y=abs(del_y-yc[i])
          id=where(del2_x lt 1 and del2_y lt 1,nidsize)
          print,nidsize,size_img
;          if nidsize gt 2 then read,idum
          if nidsize gt 1 then begin
             xmult=xc_img[id]
             ymult=yc_img[id]
             mult=1
             num=num+1
          endif
       endif
            
       if mult eq 1 then begin
          for k=0,(size(xmult))[1]-1 do begin
            
             calculate_ellipticity,xmult[k],ymult[k],es1[i],es2[i],e1,e2,infile=infile
             if (xmult[k] le nx and xmult[k] ge 0 and ymult[k] le nx and $
                 ymult[k] ge 0) then begin
                printf,1,format='(6(f10.5,1x))',xmult[k]/double(nx),$
                       ymult[k]/double(nx),e1,e2,$
                       cluster.kappa(xmult[k],ymult[k])
             endif
          endfor
       endif
       
          ;;;weak lensing
       x_pos_sing=nint(mean(idx mod nx))
       y_pos_sing=nint(mean(idx/nx))
       
       
       calculate_ellipticity,x_pos_sing,y_pos_sing,es1[i],es2[i],e1,e2,infile=infile

       if (x_pos_sing ge 0 and x_pos_sing le nx and y_pos_sing ge 0 and $
           y_pos_sing le nx) then begin
          printf,1,format='(6(f10.5,1x))',x_pos_sing/double(nx),$
                 y_pos_sing/double(nx),e1,e2,$
                 cluster.kappa(x_pos_sing,y_pos_sing)
       endif
    endif
endfor
  
close,1
close,2
print,systime()
read,idum
end

