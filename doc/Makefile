all:
	@echo Copying into HTML tree !clean manually!
	mkdir html || true
	svn export md html --force
	bash build.sh
	
publish:
	@echo "Publishing to trajectorycomputing.com"
	rsync -e "ssh -p1234" -avz html/* martin@martinwerner.de:/var/www/trajectorycomputing/doc
