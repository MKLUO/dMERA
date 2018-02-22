CC			= g++
CFLAGS		= -std=c++11 -MMD -MP
LIBS		= -luni10 
INCLUDES	= -I/usr/local/uni10/include
DBG			= -g -pg
SOURCES		= Dmera.cpp
OBJECTS		= $(SOURCES:.cpp=.o)
DEPS		= $(SOURCES:.cpp=.d) 

all: main test $(SOURCES)

main: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) main.cpp $(LIBS) $(INCLUDES) -o main.exe $(DBG)
	
test: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) test.cpp $(LIBS) $(INCLUDES) -o test.exe $(DBG)

%.o: %.cpp %.h
	$(CC) $(CFLAGS) $< $(LIBS) $(INCLUDES) -c -o $@ $(DBG)

clean:
	rm *.o *.d

-include $(DEPS)
