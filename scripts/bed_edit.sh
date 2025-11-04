#take track id from bed file and make it a column for plotting in R 

awk '
BEGIN { track_line = "" }
/^track/ { track_line = $0; next }
/^browser/ { next }
{ print $0 "\t" track_line }
' acerroh.roh.bed > acerroh.track.bed
