# Build on
# Windows 8.1 64-bit Edition
# + Cygwin 64 (2.1.0)
#  - GNU make 4.1
#  - g++ 4.9.3
#  - Boost C++ Libraries 1.58.0
#    (boost::multiprecision 1.57.0 may cause problems with Cygwin)
# + MinGW-w64 for 64 bit Windows
#  - GNU make 4.1 (mingw32-make)
#  - g++ 4.9.2
#  - Boost C++ Libraries 1.58.0
#  - Cygwin rm

TARGET = collectCards
SOURCE = collectCards.cpp
HEADER = collectCards.hpp

CXX= g++
CPPFLAGS= -std=c++11 -O2 -Wall

.PHONY: clean

all: $(TARGET)

$(TARGET): $(SOURCE) $(HEADER) Makefile
	$(CXX) $(CPPFLAGS) -o $@ $<

clean:
	rm -f $(TARGET)
