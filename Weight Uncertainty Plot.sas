/***************************************************************************************************************************
WEIGHT UNCERTAINTY PLOT COMPARING FULL SAMPLE & REPEATED HOLDOUT CHEMICAL WEIGHTS FROM WEIGHTED QUANTILE SUM REGRESSION
***************************************************************************************************************************/

*Define template;
proc template;
	define statgraph weights;
	dynamic TITLE 
			/*List of columns needed*/
			CHEMICAL_DISTRIBUTION /*All chemicals named M times (#rows = M x #chemicals)*/
			WEIGHT_DISTRIBUTION /*All M weights for each chemical (#rows = M x #chemicals)*/
			CHEMICAL /*Chemicals named for the final results (#rows = #chemicals); Note this must be in a new row*/
			WEIGHTFULL /*Final weights from full sample results (#rows = #chemicals)*/
			WEIGHTREP /*Mean weights from prepeated holdout (#rows = #chemicals)*/
			PERCENT_BAD /*#Times mean weight from prepeated holdout exceeded threshold (#rows = #chemicals)*/ 
			CHEMTYPE /*Chemical Class (#rows = #chemicals)*/
			CHEMTYPEDIST /*Chemical Class (#rows = M x #chemicals)*/;
	begingraph / pad=5 attrpriority=color;
		entrytitle textattrs=(size=9) halign=left TITLE; 
		entryfootnote textattrs=(size=8) halign=left "Notes: Bars correspond to right axis. Boxplots, diamonds, and datapoints correspond to left axis. Box plots show 25th, 50th, and 75th percentiles; whiskers show 10th and 90th percentiles. Diamonds show mean weights. Dots show all weights from 100 repeated holdouts.";
		layout overlay / xaxisopts=(display=(label tickvalues) label="Chemical" labelattrs=(Size=9 weight=normal) tickvalueattrs=(Size=8) discreteopts=(sortorder=data))
 						 yaxisopts=(label="			Weight (%)" labelattrs=(Size=9 weight=normal) 
									linearopts=(viewmin=0 viewmax=0.4 includeranges=(0-0.4) tickvaluelist=(0 0.1 0.2 0.3 0.4) tickdisplaylist=('0.0' '0.1' '0.2' '0.3' '0.4')) tickvalueattrs=(Size=9))
						 y2axisopts=(label="		# Repeated Holdouts Above Chemical of Concern Threshold" labelattrs=(Size=9 weight=normal) 
									linearopts=(viewmin=0 viewmax=100 includeranges=(0-100) tickvaluelist=(0 20 40 60 80 100) tickdisplaylist=('0' '20' '40' '60' '80' '100')) tickvalueattrs=(Size=9));
			*These are plotted underneath only in order to have correct legend;
					*Plot Mean Weight of 100 Reps;
					scatterplot y=WEIGHTREP 
								x=CHEMICAL / 
									yaxis=y 
									name="valid"
									markerattrs=(color=black size=8pt symbol=diamondfilled);
					*Plot Weight for Full Dataset;
					scatterplot y=WEIGHTFULL 
								x=CHEMICAL / 
									yaxis=y 
									name="full"
									markerattrs=(color=black size=8pt symbol=diamond);
			*Bar Chart for %Times Bad Actors over 100 Reps;
			barchartparm y=PERCENT_BAD 
						 x=CHEMICAL / 
						 	yaxis=y2 
							name="percent"
							group=CHEMTYPE
							groupdisplay=cluster
							display=(fill) 
							fillattrs=(transparency=0.9)
							baselineattrs=(thickness=0);
			*Box Plot for Distribution of Weights over 100 Reps;
			boxplot y=WEIGHT_DISTRIBUTION 
					x=CHEMICAL_DISTRIBUTION / 
							yaxis=y 
							name="box"
							group=CHEMTYPEDIST
							groupdisplay=cluster
							display=(fill median) 
							datatransparency=0.3 
							boxwidth=0.9
							fillattrs=(transparency=0.8)
							outlineattrs=(thickness=0) 
							medianattrs=(thickness=1)
							whiskerattrs=(thickness=1) 
							whiskerpercentile=10 
							capscale=0.5; 
			*Plot Weights for Each of 100 Reps;
			scatterplot y=WEIGHT_DISTRIBUTION 
						x=CHEMICAL_DISTRIBUTION / 
							yaxis=y 
							name="100 weights"
							group=CHEMTYPEDIST
							datatransparency=0.8
							markerattrs=(size=3pt symbol=circlefilled)  
							jitter=auto;
			*Plot Mean Weight of 100 Reps;
			scatterplot y=WEIGHTREP 
						x=CHEMICAL / 
							yaxis=y 
							group=CHEMTYPE
							markerattrs=(size=8pt symbol=diamondfilled);
			*Plot Weight for Full Dataset;
			scatterplot y=WEIGHTFULL 
						x=CHEMICAL / 
							yaxis=y 
							group=CHEMTYPE
							markerattrs=(size=8pt symbol=diamond);
			*Chemical of Concern threshold;
			referenceline y=0.038 / /*100/#chemicals*/
							yaxis=y 
							lineattrs=(color=black)
							curvelabel="Threshold" 
							curvelabelattrs=(Size=8 color=black) 
							curvelabelposition=min;
			*Symbol legend for mean weights;
			discretelegend "full" "valid" / 
							title="Mean Weight:" 
							titleattrs=(weight=bold size=9)
							valueattrs=(size=9)
							autoalign=(topright) 
							location=inside
							across=1;
			*Color legend for chemical groups;
			discretelegend "box" / 
							title="Chemical Type:" 
							titleattrs=(weight=normal size=9)
							valueattrs=(size=9);
		endlayout;
	endgraph;
	end;
	define style noheaderborder;
    	parent = styles.default;
		*Define colors for chemical groups so SAS doesn't assign ugly ones;
		class Graphdata1 / color=DEPPK contrastcolor=DEPPK;
		class Graphdata2 / color=BIPB contrastcolor=BIPB;
		class Graphdata3 / color=blue contrastcolor=blue;
		class Graphdata4 / color=VIBG contrastcolor=VIBG;
		class Graphdata5 / color=green contrastcolor=green;
		*Eliminate excessive borders by making them white or turning off;
    	class graphborderlines / contrastcolor=white;	
        class graphbackground / color=white ;
		style graphwalls from graphwalls / frameborder=off;
	end;
	run;

*Apply Template
ods graphics on / width=9in height=6.5in attrpriority=none border=off;
ods html style=noheaderborder;
proc sgrender data=w template=weights;
	dynamic TITLE="Chemicals of Concern Identification & Uncertainty using Traditional Full Sample vs. Repeated Holdout Validation"
			CHEMTYPEDIST="Typedist" CHEMICAL_DISTRIBUTION="Xdist" WEIGHT_DISTRIBUTION="mean_weight"  
			CHEMTYPE="Type" CHEMICAL="X" WEIGHTFULL="weight_full" WEIGHTREP="weight_rep" PERCENT_BAD="times_bad";
	run;
