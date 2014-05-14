SHELL = /bin/sh

PKGOFFSET = statlib

include $(JAVART_ROOT)/java.rules

JC = javac
JFLAGS := $(JFLAGS) -source 1.4

SRC := $(shell echo *.java)

CLASSES := $(SRC:%.java=$(CLASSESDIR)/$(PKGOFFSET)/%.class)

all:	$(CLASSES)

clean:
	rm -f $(CLASSESDIR)/$(PKGOFFSET)/*.class *.class

realclean:
	$(MAKE) clean
	rm -f *.ps *.dot *.out *~ temp*
