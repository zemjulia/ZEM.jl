function myreadlines(filename::String)
	f = open(filename)
	l = myreadlines(f)
	close(f)
	return l
end
function myreadlines(stream::IO)
	if VERSION >= v"0.6.0"
		return readlines(stream; chomp=false)
	else
		return readlines(stream)
	end
end