SHELL=/bin/sh

.SUFFIXES: .C

CPMDIR=../../opt/CPMlib-1.0.10

SRCS = \
          user_TPdomain.C \
          main.C \
        

OBJS = $(SRCS:.C=.o)

CXX      = mpic++
CXXFLAGS = `$(CPMDIR)/bin/cpm-config --cflags`
LDFLAGS  = `$(CPMDIR)/bin/cpm-config --libs`

TARGET = test


$(TARGET):		$(OBJS)
	$(CXX) $(CXXFLAGS) -o $(@) $(OBJS) $(LDFLAGS) $(LIBS)

.C.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

