function myreadlines(filename::String)
	f = open(filename)
	l = myreadlines(f)
	close(f)
	return l
end
function myreadlines(stream::IO)
	return readlines(stream; chomp=false)
end