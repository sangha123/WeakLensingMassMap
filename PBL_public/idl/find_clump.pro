pro find_clump,nx,i_pos,j_pos,xc_img,yc_img,c
  
;;; This code is used by find_images to determine the numebr of images
;;; in a group of points.
 
beta1=dblarr(nx,nx)*0.

for l=0,(size(i_pos))[1]-1 do begin
    beta1[i_pos[l],j_pos[l]]=1
endfor


x_pos=i_pos[0]
y_pos=j_pos[0]
c=0
w=1
beta_new=beta1
while (w ne 0) do begin
   index=search2d(beta_new,x_pos,y_pos,1,1,/diagonal)
   beta_new[index]=0.
   if c eq 0 then begin 
      xc_img=nint(mean(index mod nx))
      yc_img=nint(mean(index/nx))
   endif else begin 
      xc_img=[xc_img,nint(mean(index mod nx))]
      yc_img=[yc_img,nint(mean(index/nx))]
   endelse
   c=c+1
   id=where(beta_new gt 0,w)
   i0=(id mod nx)
   j0=(id/nx)
   x_pos=i0[0]
   y_pos=j0[0]
endwhile
end
