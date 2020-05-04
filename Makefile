# Makefile to compile cif2hkl.
# just type: make

# simple one-shot compile
all:	
	${FC} -O2 -o cif2hkl cif2hkl.F90 -lm
	rm *.mod

clean: 
	rm -f *.o *.mod cif2hkl

install:
	install -D cif2hkl \
		$(DESTDIR)$(prefix)/usr/bin/cif2hkl

distclean: clean

uninstall:
	-rm -f $(DESTDIR)$(prefix)/usr/bin/cif2hkl

.PHONY: all install clean distclean uninstall

