function ApplyTaper(IN,nt,nx1,tapti,taptf,tapx1i,tapx1f)

      OUT2 = zeros(size(IN));


      tt=1; tx1=1; tx2=1; tx3=1; tx4=1;

      for ix1=1:nx1

         if ix1>=1 && ix1<=tapx1i
# %            tx1 = 1 - cos(((ix1-1)/tapx1i)*pi/2);
            tx1=1/(tapx1i-1)*(ix1-1);
         end

         if ix1>=tapx1i && ix1<=nx1-tapx1f
             tx1=1;
         end

         if ix1>nx1-tapx1f && ix1<=nx1
# %           tx1 = cos(((ix1-1-nx1+tapx1f)/tapx1f)*pi/2);
		    tx1=1-1 ./(tapx1f-1)*(ix1-1-nx1+tapx1f);
         end

         for it=1:nt
               if (it>=1 && it<=tapti)
                  tt = 1.0/(tapti-1)*(it-1);
               end

               if (it>tapti && it<=nt-taptf)
                  tt = 1;
               end

               if (it>nt-taptf && it<=nt)
                  tt = 1-1.0/(taptf-1)*(it-1-nt+taptf);
		       end
		         OUT2[it,ix1] =  tt*tx1*IN[it,ix1];
         end
    end
    return OUT2
end
