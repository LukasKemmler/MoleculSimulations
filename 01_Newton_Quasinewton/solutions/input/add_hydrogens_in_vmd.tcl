package require psfgen
resetpsf
topology /Applications/VMD\ 1.9.2.app/Contents/vmd/plugins/noarch/tcl/trunctraj1.5/toppar/top_all27_lipid.rtf

#set instem 3h2o
#set instem 4h2o
#set instem 4h2o.square
set instem 3h2o.flat_triangle

segment W {
	pdb input/$instem.O.pdb
	first none
	last none
}
coordpdb input/$instem.O.pdb W

regenerate angles dihedrals
guesscoord
writepdb input/$instem.pdb
writepsf input/$instem.psf
