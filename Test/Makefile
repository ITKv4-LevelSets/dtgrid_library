DIRS = DTGrid_Test

all:
	@for i in $(DIRS); do \
	(cd $$i; $(MAKE) all); done

install:
	@for i in $(DIRS); do \
	(cd $$i; $(MAKE) install); done

clean:
	@for i in $(DIRS); do \
	(cd $$i; $(MAKE) clean); done
