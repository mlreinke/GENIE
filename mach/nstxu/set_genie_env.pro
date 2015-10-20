
; sets default environmental variables for NSTX-U

PRO set_genie_env,MACHINE=MACHINE, GENIE_PATH=GENIE_PATH, $
		  IMPSPEC_MDS_TREE=IMPSPEC_MDS_TREE, IMPSPEC_MDS_PATH=IMPSPEC_MDS_PATH, $
                  VUV_LINE_LIST=VUV_LINE_LIST, VIS_LINE_LIST=VIS_LINE_LIST

	MACHINE=getenv('MACHINE')
	if MACHINE eq '' then $
	   MACHINE='nstxu'
	
	GENIE_PATH=getenv('GENIE_PATH')
	if GENIE_PATH eq '' then $
	   GENIE_PATH='/u/mreinke/GENIE/'

	IMPSPEC_MDS_PATH=getenv('IMPSPEC_MDS_PATH')
	if IMPSPEC_MDS_PATH eq '' then $
	   IMPSPEC_MDS_PATH='\PASSIVESPEC::TOP.IMPSPEC'
	
	IMPSPEC_MDS_TREE=getenv('IMPSPEC_MDS_TREE')
	if IMPSPEC_MDS_TREE eq '' then $
	   IMPSPEC_MDS_TREE='passivespec'

	VIS_LINE_LIST=getenv('VIS_LINE_LIST')
	if VIS_LINE_LIST eq '' then $
	   VIS_LINE_LIST=GENIE_PATH+'IMPSPEC/vis_line_list.tsv'
	
	VUV_LINE_LIST=getenv('VUV_LINE_LIST')
	if VUV_LINE_LIST eq '' then $
	   VUV_LINE_LIST=GENIE_PATH+'IMPSPEC/vuv_line_list.tsv'
END

