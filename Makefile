# Makefile to compile cif2hkl.
# required: gfortran
#
# just type: make

# simple one-shot compile
all:	
	${FC} -ffree-line-length-512 \
	CFML/CFML_GlobalDeps_Linux.f90 CFML/CFML_Math_Gen.f90 \
	CFML/CFML_String_Util_gf.f90   CFML/CFML_Math_3D.f90 \
	CFML/CFML_Sym_Table.f90        CFML/CFML_Chem_Scatt.f90 \
	CFML/CFML_Symmetry.f90         CFML/CFML_Cryst_Types.f90 \
	CFML/CFML_Reflct_Util.f90      CFML/CFML_Atom_Mod.f90 \
	CFML/CFML_Geom_Calc.f90        CFML/CFML_Molecules.f90 \
	CFML/CFML_Diffpatt.f90         CFML/CFML_Magnetic_Groups.f90 \
	CFML/CFML_EisPack.f90          CFML/CFML_Form_CIF.f90 \
	CFML/CFML_Sfac.f90 -o cif2hkl cif2hkl.F90 -lm
	rm *.mod

clean: 
	rm -f *.o *.mod cif2hkl

install:
	install -D cif2hkl \
		$(DESTDIR)$(prefix)/usr/bin/cif2hkl

distclean: clean

uninstall:
	-rm -f $(DESTDIR)$(prefix)/usr/bin/cif2hkl

test:
	./cif2hkl --version
	./cif2hkl --verbose --xtal --no-outout-files examples/CaF2.cfl
	@echo Test OK

.PHONY: all install clean distclean uninstall test

