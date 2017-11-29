CC			= g++
CFLAGS		= -std=c++11
LIBS		= -luni10 -lgsl
INCLUDES	= -I/usr/local/uni10/include
DBG			= -g
SOURCES		= Dmera.cpp
OBJECTS		= $(SOURCES:.cpp=.o)

all: main $(SOURCES)

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) main.cpp $(LIBS) $(INCLUDES) -o main.exe $(DBG)
	
.cpp.o:
	$(CC) $(CFLAGS) $< $(LIBS) $(INCLUDES) -c -o $@ $(DBG)
