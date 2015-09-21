all:
	@echo "Automatic Building of Bindings is disabled in this distribution"
	@echo "Use make build to autobuild octave / python bindings"
	@echo "Alternatively:"
	@echo "make python:          Build python bindings"
	@echo "make octave:          Build octave bindings"
	@echo "sudo make python-install: (Build and) install python bindings"
	@echo "sudo make octave-install: (Build and) install octave bindings"
	@echo "sudo make install:         Install C++ header (trajcomp.hpp etc.)"
	@echo "sudo make uninstall:       Remove C++ header"
	@echo "sudo make clean:           Cleanup Development tree"
build: config-data.mk
	@test $(OCTAVE_ENABLED) && echo "Building Octave Bindings"; make -C octave || echo "Skipping Octave Bindings"
	@test $(PYTHON_ENABLED) && echo "Building Python Bindings"; make -C python || echo "Skipping Python Bindings"

config-data.mk:
	echo Trying config
	bash config.sh
python:
	make -C python
octave:
	make -C octave
python-install:
	make -C python install
octave-install:
	make -C octave install
install:
	test -d /usr/include/trajcomp || mkdir /usr/include/trajcomp
	cp doc/md/trajcomp/trajcomp.hpp /usr/include/trajcomp
	cp doc/md/trajcomp/trajcomp_files.hpp /usr/include/trajcomp
	cp rapidxml.hpp /usr/include/trajcomp
	cp rapidxml_utils.hpp /usr/include/trajcomp
	cp doc/md/trajcomp/trajcomp_geolife.hpp /usr/include/trajcomp
	cp doc/md/trajcomp/trajcomp_frechet.hpp /usr/include/trajcomp
	cp doc/md/trajcomp/trajcomp_geo.hpp /usr/include/trajcomp
	cp doc/md/trajcomp/trajcomp_all.hpp /usr/include/trajcomp
	cp doc/md/trajcomp/trajcomp_traclus.hpp /usr/include/trajcomp
	cp bloom/murmur.hpp /usr/include/trajcomp
	cp bloom/bloomfilter.hpp /usr/include/trajcomp
	cp bloom/countingbloomfilter.hpp /usr/include/trajcomp
	cp bloom/timedecayingbloomfilter.hpp /usr/include/trajcomp
uninstall:
	rm -f /usr/include/trajcomp/trajcomp.hpp
	rm -f /usr/include/trajcomp/trajcomp_files.hpp
	rm -f /usr/include/trajcomp/rapidxml.hpp
	rm -f /usr/include/trajcomp/rapidxml_utils.hpp
	rm -f /usr/include/trajcomp/trajcomp_geolife.hpp
	rm -f /usr/include/trajcomp/trajcomp_frechet.hpp
	rm -f /usr/include/trajcomp/trajcomp_geo.hpp
	rm -f /usr/include/trajcomp/trajcomp_all.hpp
	rm -f /usr/include/trajcomp/trajcomp_traclus.hpp	
	rm -f /usr/include/trajcomp/murmur.hpp
	rm -f /usr/include/trajcomp/bloomfilter.hpp
	rm -f /usr/include/trajcomp/countingbloomfilter.hpp
	rm -f /usr/include/trajcomp/timedecayingbloomfilter.hpp
	rmdir /usr/include/trajcomp
clean:
	make -C octave clean
	make -C python clean
	rm -f config-data.mk
