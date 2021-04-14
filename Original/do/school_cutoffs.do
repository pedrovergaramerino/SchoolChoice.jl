* SCHOOL CUTOFFS

capture program drop Cutoffs_graph
	program  Cutoffs_graph
	syntax [, infile(string) ceiling(string) step(string) outfile(string) samp(string) s2(string)]

	   gen empty=0		
	   
	   replace School1=. if School1==0
	   replace School2=. if School2==0
	   replace School3=. if School3==0
	   replace School4=. if School4==0
	   replace School5=. if School5==0
	   replace School6=. if School6==0				
								
	   if "`samp'"=="100"{
			local y1=9
			local x1=0.06
			local y2=12
			local x2=0.14
			local y3=8
			local x3=0.23			
			local y4=8
			local x4=0.42			
			local y5=7
			local x5=0.58		  
			local y6=6
			local x6=0.77
			}				
	   if "`samp'"=="200"{
			local y1=12
			local x1=0.06
			local y2=15
			local x2=0.14
			local y3=11
			local x3=0.23			
			local y4=8
			local x4=0.42			
			local y5=7
			local x5=0.58		  
			local y6=6
			local x6=0.77
			}				
	   if "`samp'"=="500"{
			local y1=18
			local x1=0.06
			local y2=22
			local x2=0.14
			local y3=15
			local x3=0.23			
			local y4=13
			local x4=0.42			
			local y5=9
			local x5=0.58		  
			local y6=7.5
			local x6=0.75
			}		
	   if "`samp'"=="5000"{
			local y1=45.5
			local x1=0.08
			local y2=42
			local x2=0.18
			local y3=38
			local x3=0.25			
			local y4=37
			local x4=0.42			
			local y5=26.5
			local x5=0.55		  
			local y6=25
			local x6=0.7
			}			
	   twoway (line School1 Cutoff, lw(medthick) lp(solid) lcolor(navy)) ///	   
	   (line School5 Cutoff, lw(medthick) lp(solid) lcolor(purple)) ///
	   (line School2 Cutoff,  lw(medthick) lp(solid) lcolor(green)) ///	 
  	   (line School3 Cutoff,  lw(medthick) lp(solid) lcolor(maroon)) ///
	   (line School4 Cutoff,  lw(medthick) lp(solid) lcolor(ebblue)) ///	   
	   (line School6 Cutoff, lw(medthick) lp(solid) lcolor(dkorange)) ///
	   (line empty Cutoff,  lw(thin) lp(solid) lcolor(black)) ///
	   (scatteri `s2', recast(line) lw(vthin) color(black) legend(off)), ///
	   xtitle("School cutoff", size(medlarge) height(5)) ///
	   ytitle("Density", size(medlarge) height(5)) ///
	   xscale(range(0 1)) ///
	   xlabel(0(0.1)1, format(%2.1fc)) ///
	   yscale(range(0 `ceiling')) ///
	   ylabel(0(`step')`ceiling', glcolor(gs14) glwidth(thin) gmin gmax) ///
	   legend(off) ///
	   graphregion( fcolor(white) lcolor(white)) ///	
	   yline(0, lw(thin) lcolor(black)) ///
	   text(`y1' `x1' "School 1", color(navy)) ///
	   text(`y2' `x2' "School 2", color(green)) ///
	   text(`y3' `x3' "School 5", color(purple)) ///
	   text(`y4' `x4' "School 6", color(dkorange)) ///
	   text(`y5' `x5' "School 4", color(ebblue)) ///
	   text(`y6' `x6' "School 3", color(maroon))   

end

capture program drop Cutoffs_graph_new
	program  Cutoffs_graph_new
	syntax [, ceiling(string) step(string) samp(string) s2(string)]

	   gen empty=0		
	   
	   replace School1=. if School1==0
	   replace School2=. if School2==0
	   replace School3=. if School3==0
	   replace School4=. if School4==0
	   replace School5=. if School5==0
	   replace School6=. if School6==0				

	   	if "`samp'"=="500"{
			local y1=17
			local x1=0.08
			local y2=15.8
			local x2=0.18
			local y3=14.3
			local x3=0.23			
			local y4=11.7
			local x4=0.423			
			local y5=8
			local x5=0.58		  
			local y6=7.3
			local x6=0.72
			}
			
	   if "`samp'"=="100"{
			local y1=7.9
			local x1=0.06
			local y2=7
			local x2=0.18
			local y3=6
			local x3=0.27			
			local y4=6.2
			local x4=0.45			
			local y5=4.4
			local x5=0.58		  
			local y6=3.7
			local x6=0.79
			}
	   twoway (line School1 Cutoff, lw(medthick) lp(solid) lcolor(navy)) ///	   
	   (line School5 Cutoff, lw(medthick) lp(solid) lcolor(purple)) ///
	   (line School2 Cutoff,  lw(medthick) lp(solid) lcolor(green)) ///	 
  	   (line School3 Cutoff,  lw(medthick) lp(solid) lcolor(maroon)) ///
	   (line School4 Cutoff,  lw(medthick) lp(solid) lcolor(ebblue)) ///	   
	   (line School6 Cutoff, lw(medthick) lp(solid) lcolor(dkorange)) ///
	   (line empty Cutoff,  lw(thin) lp(solid) lcolor(black)) ///
	   (scatteri `s2', recast(line) lw(vthin) color(black) legend(off)), ///
	   xtitle("School cutoff", size(medlarge) height(5)) ///
	   ytitle("Density", size(medlarge) height(5)) ///
	   xscale(range(0 1)) ///
	   xlabel(0(0.1)1, format(%2.1fc)) ///
	   yscale(range(0 `ceiling')) ///
	   ylabel(0(`step')`ceiling', glcolor(gs14) glwidth(thin) gmin gmax) ///
	   legend(off) ///
	   graphregion( fcolor(white) lcolor(white)) ///	
	   yline(0, lw(thin) lcolor(black)) ///
	   text(`y1' `x1' "School 1", color(navy)) ///
	   text(`y2' `x2' "School 2", color(green)) ///
	   text(`y3' `x3' "School 5", color(purple)) ///
	   text(`y4' `x4' "School 6", color(dkorange)) ///
	   text(`y5' `x5' "School 4", color(ebblue)) ///
	   text(`y6' `x6' "School 3", color(maroon))   

end

