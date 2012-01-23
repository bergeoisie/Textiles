OBJS = textileHelper.o
CC = g++
CFLAGS = -c -O3
LFLAGS = -O3
INCS = -I ./boost_1_47_0/ -I ./TNT -I ./JAMA
gmcube: $(OBJS) 
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o gmcube gmcube.cpp

gmsquare: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o gmsquare gmsquare.cpp

gmctwo: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o gmctwo gmctwo.cpp
	
gmNew21: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o gmNew21 gmNew21.cpp
	
allNew: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o gmNew21 gmNew21.cpp
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o gmNew231 gmNew231.cpp
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o gmNew2341 gmNew2341.cpp
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o gmNew2143 gmNew2143.cpp	
 
junk: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) $(INCS) -o junk junk.cpp

textileHelper.o: textileHelper.cpp
	$(CC) $(CFLAGS) $(INCS) textileHelper.cpp

clean: 
	\rm *.o
	\rm gmcube
	\rm gmsquare
	\rm gmctwo
	\rm junk