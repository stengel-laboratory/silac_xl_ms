# silac_xl_ms
Files for evaluating the SILAC XL-MS workflow with the xQuest software.

This repository contains the modified xQuest (source) files for enabling SILAC XL-MS analysis.
See: ref. to paper once published.

Get xQuest here: http://proteomics.ethz.ch/cgi-bin/xquest2_cgi/download.cgi
xQuest references: doi:10.1038/nmeth.1192; doi:10.1038/nmeth.2103; doi:10.1038/nprot.2013.168

Instructions:

1) Copy the provided xQuest folder structure including the files over your installation (or rather make a 2nd installation)
2) Begin with the regular xQuest analysis pipeline


Details:
The following files were modified:

- xquest/deffiles/xQuest/modules/Xquest_Digest.pm: 
	- Added a customized trypsin variant (enzyme #99) allowing to cut the B,J and U amino acids
- xquest/deffiles/mass_table.def:
	- Added the three unnatural amino acids B,J and U
- xquest/deffiles/xQuest/xquest.def:
	- Modified to use our customized trypsin enzyme

Note that our modification is based on xQuest version 2.1.3
