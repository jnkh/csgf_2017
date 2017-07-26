all:
	pgc++ -g -acc -O2 mb.cc -o mb

gpu:
	pgc++ -ta=tesla:cc60 -Minfo=accel -acc mb.cc -o mb-gpu

clean:
	${RM} mb
