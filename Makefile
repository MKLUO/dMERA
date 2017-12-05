CC			= g++
CFLAGS		= -std=c++11
LIBS		= -luni10 
INCLUDES	= -I/usr/local/uni10/include
DBG			= -g
SOURCES		= Dmera.cpp
OBJECTS		= $(SOURCES:.cpp=.o)

all: main $(SOURCES) test

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) main.cpp $(LIBS) $(INCLUDES) -o main.exe $(DBG)
	
test: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test.cpp $(LIBS) $(INCLUDES) -o test.exe $(DBG)

.cpp.o:
	$(CC) $(CFLAGS) $< $(LIBS) $(INCLUDES) -c -o $@ $(DBG)

clear:
	rm *.o
