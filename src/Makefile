all:
	pgc++ -g -O2 mb.cc -o mb

gpu:
	pgc++ -O2 -ta=tesla:cc60 -Minfo=accel -acc mb.cc -o mb-gpu

clean:
	${RM} mb
