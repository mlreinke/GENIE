; sets default environmental variables for C-Mod

PRO get_genie_env

	MACHINE=getenv('MACHINE')
	if MACHINE eq '' then $
	   MACHINE='cmod'
	
	GENIE_PATH=getenv('GENIE_PATH')
	if GENIE_PATH eq '' then $
	   GENIE_PATH='/usr/local/cmod/idl/GENIE/'

	IMPSPEC_MDS_PATH=getenv('IMPSPEC_MDS_PATH')
	if IMPSPEC_MDS_PATH eq '' then $
	   IMPSPEC_MDS_PATH='\SPECTROSCOPY::TOP.IMPSPEC'
	
	IMPSPEC_MDS_TREE=getenv('IMPSPEC_MDS_TREE')
	if IMPSPEC_MDS_TREE eq '' then $
	   IMPSPEC_MDS_TREE='spectroscopy'

	VIS_LINE_LIST=getenv('VIS_LINE_LIST')
	if VIS_LINE_LIST eq '' then $
	   VIS_LINE_LIST=GENIE_PATH+'IMPSPEC/vis_line_list.tsv'
	
	VUV_LINE_LIST=getenv('VUV_LINE_LIST')
	if VUV_LINE_LIST eq '' then $
	   VUV_LINE_LIST=GENIE_PATH+'IMPSPEC/vuv_line_list.tsv'
END

