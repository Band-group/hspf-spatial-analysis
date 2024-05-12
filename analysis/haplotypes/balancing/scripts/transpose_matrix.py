import sys
rows = 0
cols = None
data = []
for line in sys.stdin:
	a = line.strip().split( "\t" )
	if cols is None:
		cols = len(a)
	elif len(a) != cols:
		raise Exception( "Line %d has %d cols, expected %d." % ( rows+1, len(a), cols ))
	data += a
	rows += 1

#print( "%d x %d" % ( rows, cols ))

for i in range(0,cols):
	start = i
	end = i+rows*cols
	by = cols
	sys.stdout.write( ' '.join( data[start:end:by] ) + "\n" )
