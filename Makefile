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