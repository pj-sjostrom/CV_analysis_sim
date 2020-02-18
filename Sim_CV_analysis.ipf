#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Menu "Macros"
	"Set up simulations",SetUpSimulations()
End

// Set up the simulations

Function SetUpSimulations()

	Variable/G	simLTP = 0		// Boolean: Simulate LTP (otherwise LTD)
	String/G		caseStr = ""	// List of possible choices of simulations
	caseStr += "1 - pre expression with increasing noise levels;"
	caseStr += "2 - post expression with increasing noise levels;"
	caseStr += "3 - increasing outlier, pre expression;"
	caseStr += "4 - increasing outlier, post expression;"
	caseStr += "5 - increasing baseline trend, pre expression;"
	caseStr += "6 - increasing baseline trend, post expression;"
	caseStr += "7 - graded development of pre expression;"
	caseStr += "8 - graded development of post expression;"
	String/G		MostRecentSimulation = "<null>"		// Remember what the most recent simulation was (text)
	Variable/G	LastSimNum = 3		// Remember what the most recent simulation was (number)
	Variable/G	p_before = NaN		// Probability of release before induction of plastcity
	Variable/G	p_after = NaN		// Probability of release after
	Variable/G	q_before = NaN		// Quantal amplitude in mV terms before
	Variable/G	q_after = NaN		// Quantal amplitude in mV terms after
		
	DoWindow/K SimulationPanel		//	Redraw the simulation panel (ensure that restart is possible)
	Execute "SimulationPanel()"

end

// Draw the simulation panel

Window SimulationPanel() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(4,53,412,168)
	PopupMenu simModePopup,pos={4.00,1.00},size={400-250,23.00},bodyWidth=400,proc=PopExecProc,title="Pick a simulation"
	PopupMenu simModePopup,mode=0,value= #"caseStr",fSize=10
	CheckBox LTPcheck,pos={4.00,24.00},size={400,16.00},title="LTP (otherwise LTD)"
	CheckBox LTPcheck,fSize=10,variable= simLTP,proc=LTPCheckProc
	SetVariable LastSimWas,pos={4.00,48.00},size={400,16.00},title="Last sim"
	SetVariable LastSimWas,fSize=10,value= MostRecentSimulation,noedit=1
	ValDisplay valdisp0,pos={4.00,72.00},size={196.00,15.00},title="p before"//,fColor=(65535/2,65535/2,65535/2),frame=0
	ValDisplay valdisp0,fSize=10,limits={0,0,0},barmisc={0,1000},value= #"p_before"
	ValDisplay valdisp1,pos={204.00,72.00},size={200.00,15.00},title="q before (mV)"
	ValDisplay valdisp1,fSize=10,limits={0,0,0},barmisc={0,1000},value= #"q_before"
	ValDisplay valdisp2,pos={4.00,96.00},size={196.00,15.00},title="p after"
	ValDisplay valdisp2,fSize=10,limits={0,0,0},barmisc={0,1000},value= #"p_after"
	ValDisplay valdisp3,pos={204.00,96.00},size={200.00,15.00},title="q after (mV)"
	ValDisplay valdisp3,fSize=10,limits={0,0,0},barmisc={0,1000},value= #"q_after"
EndMacro

// Execute simulation when popup is selected

Function PopExecProc(pa) : PopupMenuControl
	STRUCT	WMPopupAction &pa
	
	NVAR		LastSimNum

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			LastSimNum = popNum
			reNoise(popNum)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

// Execute simulation when LTP/LTD checkbox is altered

Function LTPCheckProc(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	NVAR		LastSimNum

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			reNoise(LastSimNum)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

// Carry out simulation numerous times for each setting, to account for noise

Function reNoise(which)
	Variable	which

	SVAR		caseStr
	NVAR		simLTP
	
	Variable	n = 150
	Variable	Noise1 = 0.2e-3
	Variable nSteps = 6
	Variable theNoise = 1e-3
	Variable	i,j
	
	Variable/G	outlierVal = 0
	Variable/G	trendVal = 0
	Variable/G	gradeVal = 0
	
	String	legStr = "empty"
	
	String	grList = ""
	
	j = 0
	do
		switch(which)			// Cases as described by caseStr above in SetUpSimulations
			case 1:
				theNoise = Noise1/2^(nSteps-j-1)
				if (j==0)
					theNoise = 0
				endif
				legStr = "Noise sdev:\r"+num2str(theNoise*1e3)+" mV"
				break
			case 2:
				theNoise = Noise1/2^(nSteps-j-1)
				if (j==0)
					theNoise = 0
				endif
				legStr = "Noise sdev:\r"+num2str(theNoise*1e3)+" mV"
				break
			case 3:
				theNoise = 0.1e-3
				outlierVal = 0.1e-3*2^j
				if (j==0)
					outlierVal = 0
				endif
				legStr = "Outlier: "+num2str(outlierVal*1e3)+" mV"
				break
			case 4:
				theNoise = 0.1e-3
				outlierVal = 0.1e-3*2^j
				if (j==0)
					outlierVal = 0
				endif
				legStr = "Outlier:\r"+num2str(outlierVal*1e3)+" mV"
				break
			case 5:
				theNoise = 0.1e-3
				trendVal = 0.6e-6*2^j
				if (j==0)
					trendVal = 0
				endif
				legStr = "Trend:\r"+num2str(trendVal*1e6*60/10)+" µV/min"
				break
			case 6:
				theNoise = 0.1e-3
				trendVal = 0.6e-6*2^j
				if (j==0)
					trendVal = 0
				endif
				legStr = "Trend:\r"+num2str(trendVal*1e6*60/10)+" µV/min"
				break
			case 7:
				theNoise = 0.1e-3
				gradeVal = j/(nSteps-1)
				legStr = "Graded:\r"+num2str(gradeVal)+""
				break
			case 8:
				theNoise = 0.1e-3
				gradeVal = j/(nSteps-1)
				legStr = "Graded:\r"+num2str(gradeVal)+""
				break
			default:
				theNoise = 0.1e-3
		endswitch
		Make/O/N=(0) gatherCVinvSq,gatherMeanNorm
		i = 0
		do
			noiseSim(theNoise,which)
			WAVE	wCVinvSq
			WAVE	wMeanNorm
			gatherCVinvSq[numpnts(gatherCVinvSq)] = {wCVinvSq[1]}
			gatherMeanNorm[numpnts(gatherMeanNorm)] = {wMeanNorm[1]}
			i += 1
		while(i<n)
	
		makeEnsCVplot(j+1)
		legend/A=LT/F=0/B=1/X=0.00/Y=0.00 legStr
	
		grList += "ensCVplot_"+num2str(j+1)+";"

		j += 1
	
	while(j<nSteps)

	showData()
	ArrangeGraphs(";thePlot;",3,3)
	ArrangeGraphs(";;;;;;"+grList,3,6)
	MoveGraphsToFront(grList)
	
	if ((which==3) %| (which==4))
		Fix3()
	endif
	if ((which==5) %| (which==6))
		Fix5()
	endif

	SVAR	MostRecentSimulation
	if (simLTP)
		MostRecentSimulation = "LTP: "
	else
		MostRecentSimulation = "LTD: "
	endif	
	MostRecentSimulation += stringFromList(which-1,caseStr)+" at "+Time()
	print MostRecentSimulation

End

// Carry out simulation once, to be called numerous times

Function noiseSim(theNoise,which)
	Variable	theNoise
	Variable which
	
	NVAR		outlierVal
	NVAR		trendVal
	NVAR		gradeVal

	NVAR		simLTP
	
	Variable nBefore = 60
	Variable	nAfter = nBefore*4
	Variable n = 5

	// Pre
	NVAR p_before
	NVAR p_after
	NVAR q_before
	NVAR q_after
	if (simLTP)
		p_before = 0.4
		p_after = 0.55
	else
		p_before = 0.55
		p_after = 0.4
	endif
	q_before = 0.35
	q_after = 0.35
	// Post
	if ((which==2) %| (which==4) %| (which==6) %| (which==8))
		p_before = 0.5
		p_after = 0.5
		if (simLTP)
			q_before = 0.28
			q_after = 0.4
		else
			q_before = 0.4
			q_after = 0.28
		endif
	endif
	
	Make/O/N=(nBefore) wBefore
	Make/O/N=(nAfter) wAfter
	WAVE	wBefore
	WAVE	wAfter
	
	SetScale/P x 0,0.1666666,"min", wBefore,wAfter
	SetScale d 0,0,"V", wBefore,wAfter
	
	wBefore = simBin(p_before,n)
	wBefore *= q_before*1e-3
	wBefore += gNoise(theNoise)
	Duplicate/O wBefore,wBefore_Orig
	// Add baseline outlier
	wBefore[10] += outlierVal
	//	Add baseline linear trend
	wBefore += trendVal*p-trendVal*nBefore/2

	wAfter = simBin(p_after,n)
	Duplicate/O wAfter,wAfter_Orig
	if (which==7)
		wAfter = simBin(gradeVal*(p_before-p_after)*exp(-p/100)+p_after,n)
	endif
	wAfter *= q_after*1e-3
	if (which==8)
		wAfter = simBin(p_after,n)*(gradeVal*(q_before*1e-3-q_after*1e-3)*exp(-p/100)+q_after*1e-3)
	endif
	wAfter += gNoise(theNoise)
	
	doCVanalysis(theNoise)

End

// Carry out CV analysis on a simulated experiment

Function doCVanalysis(theNoise)
	Variable	theNoise

	WAVE	wBefore
	WAVE	wAfter
	
	WaveStats/Q wBefore
	Variable	m1 = V_avg
	Variable	s1 = V_SDev
	Variable CV1 = s1/m1
	
	WaveStats/Q wAfter
	Variable	m2 = V_avg
	Variable	s2 = V_SDev//-theNoise
	Variable CV2 = s2/m2
	
	Make/O/N=(2) wCVinvSq,wMeanNorm,wDiagX,wDiagY
	wDiagX = {0,4.0001}
	wDiagY = {0,4}
	wCVinvSq = 1
	wMeanNorm = 1
	wCVinvSq[1] = (CV1/CV2)^2
	wMeanNorm[1] = m2/m1

End

// Produce ensemble average CV analysis data with error bars and plot that along with CV analysis of individual experiments

Function makeEnsCVplot(j)
	variable j
	
	WAVE	wCVinvSq
	WAVE	wMeanNorm
	
	WAVE	gatherCVinvSq
	WAVE	gatherMeanNorm
	
	Duplicate/O wCVinvSq,$("wCVinvSq_"+num2str(j))
	Duplicate/O wMeanNorm,$("wMeanNorm_"+num2str(j))
	
	WAVE	wCVinvSq_j = $("wCVinvSq_"+num2str(j))
	WAVE	wMeanNorm_j = $("wMeanNorm_"+num2str(j))

	Duplicate/O gatherCVinvSq,$("gatherCVinvSq_"+num2str(j))
	Duplicate/O gatherMeanNorm,$("gatherMeanNorm_"+num2str(j))
	
	doWindow/K $("ensCVplot_"+num2str(j))
	Display $("wCVinvSq_"+num2str(j)) vs $("wMeanNorm_"+num2str(j)) as "CV analysis "+num2str(j)
	doWindow/C $("ensCVplot_"+num2str(j))
	ModifyGraph mode=4,marker=19,rgb=(14080,30464,47872)//,msize=6
	SetAxis/A/E=1 left
	SetAxis/A/E=1 bottom
	ModifyGraph width={Plan,1,bottom,left}
	ModifyGraph manTick(left)={0,1,0,0},manMinor(left)={9,5},manTick(bottom)={0,1,0,0},manMinor(bottom)={9,5}
	DoUpdate
	Label left,"1/CV\S2\M (norm)"
	Label bottom,"mean (norm)"
	
	AppendToGraph $("gatherCVinvSq_"+num2str(j)) vs $("gatherMeanNorm_"+num2str(j))
	ModifyGraph mode($("gatherCVinvSq_"+num2str(j)))=3,marker($("gatherCVinvSq_"+num2str(j)))=8,rgb($("gatherCVinvSq_"+num2str(j)))=(48059,48059,48059)
	
	WAVE	wDiagY
	WAVE	wDiagX
	AppendToGraph wDiagY vs wDiagX
	ModifyGraph lstyle(wDiagY)=11,rgb(wDiagY)=(30583,30583,30583)
	ReorderTraces $("wCVinvSq_"+num2str(j)),{$("gatherCVinvSq_"+num2str(j)),wDiagY}
	
	Make/O/N=2 $("gatherCVinvSq_SEM_"+num2str(j)),$("gatherMeanNorm_SEM_"+num2str(j))
	WAVE	gatherCVinvSq_SEM_j = $("gatherCVinvSq_SEM_"+num2str(j))
	WAVE	gatherMeanNorm_SEM_j = $("gatherMeanNorm_SEM_"+num2str(j))
	gatherCVinvSq_SEM_j = 0
	gatherMeanNorm_SEM_j = 0
	WaveStats/Q $("gatherCVinvSq_"+num2str(j))
	wCVinvSq_j[1] = V_avg
	gatherCVinvSq_SEM_j[1] = V_SDev
	
	WaveStats/Q gatherMeanNorm
	wMeanNorm_j[1] = V_avg
	gatherMeanNorm_SEM_j[1] = V_SDev
	
	ErrorBars/W=$("ensCVplot_"+num2str(j)) $("wCVinvSq_"+num2str(j)) XY,wave=(gatherMeanNorm_SEM_j,gatherMeanNorm_SEM_j),wave=(gatherCVinvSq_SEM_j,gatherCVinvSq_SEM_j)

	SetDrawLayer UserBack
	SetDrawEnv xcoord= prel,ycoord= left,dash= 1
	DrawLine 0,1,1,1
	SetDrawEnv xcoord= bottom,ycoord= prel,dash= 1
	DrawLine 1,0,1,1
	SetDrawLayer UserFront

	NVAR	simLTP
	if (simLTP)
		SetAxis left,0,4
		SetAxis bottom,0,4
	else
		SetAxis left,0,2
		SetAxis bottom,0,2
	endif
	
	doUpdate

End

// Show a single simulated LTP or LTD experiment

Function showData()

	WAVE	wBefore
	WAVE	wAfter
	
	Variable	afterShift = 20
	
	Variable	binSize = numpnts(wBefore)/2
	Variable	nBins = 2+1+2*4
	Make/O/N=(nBins) LTP_wave,time_wave,SEM_wave
	Make/O/N=(1) nanWave = NaN
	Variable	i = 0
	Variable j
	do
		if (i<2)
			WAVE	theWave = wBefore
			j = 0
		else
			if (i>2)
				WAVE	theWave = wAfter
				j = -90
			else
				WAVE	theWave = nanWave
				j = 0
			endif
		endif
		WaveStats/Q/R=[0+i*numpnts(wBefore)/2+j,numpnts(wBefore)/2-1+i*numpnts(wBefore)/2+j] theWave
		LTP_wave[i] = V_avg
		SEM_wave[i] = V_SEM
		time_wave[i] = 0.1666666*(numpnts(wBefore)/4+i*numpnts(wBefore)/2)
		if (i>2)
			time_wave[i] += 10-0.1666666*numpnts(wBefore)/2
		endif
		i += 1
	while (i<nBins)

	doWindow/K thePlot
	Display wBefore,wAfter as "The Data"
	doWindow/C thePlot
	ModifyGraph mode=3,marker=8,opaque=1,rgb=(14080,30464,47872)
	ModifyGraph offset(wAfter)={afterShift,0}
	
	Variable	lThick = 1
	AppendToGraph LTP_wave vs time_wave
	ModifyGraph mode(LTP_wave)=4,marker(LTP_wave)=19,lsize(LTP_wave)=lThick,rgb(LTP_wave)=(0,0,0)
	ErrorBars/T=(lThick)/L=(lThick) LTP_wave Y,wave=(SEM_wave,SEM_wave)
	
	SetDrawLayer userback
	SetDrawEnv xcoord= bottom,ycoord= left,dash= 1
	DrawLine 0,mean(wBefore),rightx(wAfter)+afterShift,mean(wBefore)

	SetDrawEnv xcoord= bottom,fillfgc= (56797,56797,56797),linethick= 0.00
	DrawRect 10,0,20,1
	
	SetDrawEnv xcoord= bottom,textxjust= 1,textyjust= 0,textrot= 0,fsize= 10
	DrawText 15,1,"induction"
	
End

// Simulate the binomial distribution

Function simBin(p,N)
	Variable	p,N
	
	Variable i = 0
	Variable	outcome = 0
	
	do
		if (eNoise(0.5)+0.5<p)
			outcome += 1
		endif
		i += 1
	while (i<N)
	
	Return outcome
	
End

// To organize plots, need to know the screen size

Function GetScreenSize(yAxis)
	Variable	yAxis

	String		ScreenSizeStr = IgorInfo(0)
	ScreenSizeStr = ScreenSizeStr[StrSearch(ScreenSizeStr,"RECT=",0),inf]
	
	Variable	xSize = Str2Num(StringFromList(2,ScreenSizeStr,","))
	Variable	ySize = Str2Num(StringFromList(3,ScreenSizeStr,","))
	
	if (yAxis)
		Return ySize
	else
		Return xSize
	endif

End

// Organize plots in rows and columns

Function ArrangeGraphs(ListOfGraphNames,Rows,Columns)		// Put graphs in nice rows and columns
	String		ListOfGraphNames						// Semi-colon separated list
	Variable	Rows
	Variable	Columns
	
	Variable	xSize = GetScreenSize(0)
	Variable	ySize = GetScreenSize(1)-20
	Variable	xSkip = 16
	Variable	ySkip = 40

	Variable	xPos = 8
	Variable	yPos = 64

	Variable	width = (xSize-xPos)/Columns-xSkip
	Variable	height = (ySize-yPos)/Rows-ySkip

	Variable	ScSc = 72/ScreenResolution				// Screen resolution
	
	Variable	x1
	Variable	y1
	String		currName
	
	Variable	i
	Variable	r = 0
	Variable	c = 0
	Variable	n = ItemsInList(ListOfGraphNames)
	i = 0
	do
		currName = StringFromList(i,ListOfGraphNames)
		if (!(StringMatch(currName,"")))
			x1 = xPos+c*(width+xSkip)
			y1 = yPos+r*(height+ySkip)
			DoWindow $(currName)			// Only try to move windows that exist
			if (V_flag>0)
			 	MoveWindow/W=$(currName) x1*ScSc, y1*ScSc ,(x1+width)*ScSc,(y1+height)*ScSc
			endif
		 endif
	 	c += 1
	 	if (c>Columns-1)
	 		c = 0
	 		r += 1
	 		if (r>Rows-1)
	 			r = 0
	 		endif
	 	endif
		i += 1
	while(i<n)

End

// Move graphs to front

Function MoveGraphsToFront(ListOfGraphNames)		// Graphs to front
	String		ListOfGraphNames						// Semi-colon separated list

	String		currName
	
	Variable	i
	Variable	n = ItemsInList(ListOfGraphNames)
	i = 0
	do
		currName = StringFromList(i,ListOfGraphNames)
		if (!(StringMatch(currName,"")))
			DoWindow/F $(currName)
		endif
		i += 1
	while(i<n)

End

// Case 3 (and 4) needs a special cosmetic tweak to look like it looks in the paper

Function Fix3()

	DoWindow/F thePlot
	WAVE wBefore,wBefore_Orig
	ModifyGraph sep(left)=1,manTick(left)={0,1,-3,0},manMinor(left)={1,50}	
	ModifyGraph rgb(wBefore[10])=(60652,36494,37265),marker(wBefore[10])=19,mrkStrokeRGB(wBefore[10])=(58652,0,7208)	
	
	AppendToGraph wBefore_Orig
	Variable recall = wBefore_Orig[10]
	wBefore_Orig = NaN
	wBefore_Orig[10] = recall
	ModifyGraph mode(wBefore_Orig)=3
	ModifyGraph marker(wBefore_Orig)=19
	ModifyGraph rgb(wBefore_Orig)=(60652,36494,37265),mrkStrokeRGB(wBefore_Orig[10])=(58652,0,7208)	

	SetDrawLayer/K UserFront
	SetDrawEnv xcoord= bottom,ycoord= left, arrow= 1,arrowfat= 1.00,linethick= 2.00
	DrawLine 1.666666666666,wBefore[10]-0.0032*0.9,1.666666666666,wBefore[10]-0.0032*0.1
	
	SetDrawEnv xcoord= bottom,ycoord= left,textyjust= 1,fsize= 10,textrot= 90
	DrawText 2,wBefore[10]-0.0032*0.5,"+3.2 mV"

	SetAxis/N=1 left -0.0005,*
	SetAxis/A/N=1 bottom

End

// Case 5 (and 6) needs a special cosmetic tweak to look like it looks in the paper

Function Fix5()

	DoWindow/F thePlot
	WAVE	wBefore,wBefore_Orig

	SetAxis/A/N=1 left
	SetAxis/A/N=1 bottom
	
End
